import sys
import time
import argparse

import numpy as np
import pandas as pd
import anndata as ad
import ruptures as rpt
from tqdm import tqdm


### Helper functions ###
def sparse_band_corr(X, k):
    """
    Compute only the first k diagonals of the correlation matrix of X,
    stored in sparse format.

    Parameters
    ----------
    X : scipy.sparse matrix, np.ndarray
        Input data matrix of shape (n_samples, n_features).
    k : int
        Number of diagonals to compute.

    Returns
    -------
    band_corr : 2D numpy array
        The banded correlation matrix of shape (n_features, n_features).
    """
    n, p = X.shape

    band_corr = np.zeros((p, p), dtype=np.float32)

    for block in tqdm(range(0, p // k)):
        start = block * k
        end = min(start + k, p)

        if hasattr(X, "toarray"):
            X_block = X[:, start:end].toarray()
            off_block = X[:, end : 2 * end - start].toarray()
        elif isinstance(X, np.ndarray):
            X_block = X[:, start:end]
            off_block = X[:, end : 2 * end - start]
        else:
            raise ValueError("Input X must be a scipy.sparse matrix or a numpy ndarray.")

        Xc = X_block - X_block.mean(axis=0)
        Yc = off_block - off_block.mean(axis=0)

        X_std = np.sqrt((Xc**2).sum(axis=0))
        Y_std = np.sqrt((Yc**2).sum(axis=0))
        X_std[X_std == 0] = np.inf
        Y_std[Y_std == 0] = np.inf

        block_corr = np.dot(Xc.T, Xc) / np.outer(X_std, X_std)
        off_block_corr = np.dot(Xc.T, Yc) / np.outer(X_std, Y_std)

        triu_ind = np.triu_indices(n=off_block_corr.shape[0], k=1, m=off_block_corr.shape[1])
        off_block_corr[triu_ind] = 0

        band_corr[start:end, start:end] = block_corr
        band_corr[start:end, end : 2 * end - start] = off_block_corr
        band_corr[end : 2 * end - start, start:end] = off_block_corr.T

    # Last block
    start = p - k
    end = p
    X_block = X[:, start:end].toarray()

    Xc = X_block - X_block.mean(axis=0)
    X_std = np.sqrt((Xc**2).sum(axis=0))
    X_std[X_std == 0] = np.inf
    block_corr = np.dot(Xc.T, Xc) / np.outer(X_std, X_std)
    band_corr[start:end, start:end] = block_corr

    return band_corr


def distance_kernel(sigma, truncate=4.0, radius=None):
    """
    Create a square Gaussian distance kernel.

    Parameters
    ----------
    sigma : float
        Standard deviation of the Gaussian.
    truncate : float, optional
        Truncate the kernel at this many standard deviations. Default is 4.0.
    radius : int, optional
        Radius of the kernel. If None, it is set to int(truncate * sigma).

    Returns
    -------
    kernel : 2D numpy array
        The Gaussian distance kernel.
    """
    if radius is None:
        radius = int(truncate * sigma)

    width = 1 + 2 * radius
    kernel = np.zeros((width, width))
    for offset in range(1, width):
        kernel += offset * (np.eye(width, k=offset) + np.eye(width, k=-offset))

    kernel += np.rot90(kernel)
    kernel /= 2
    kernel = np.exp(-(kernel**2) / (2 * sigma**2))
    # np.fill_diagonal(kernel, 0) # To not count self-correlation
    kernel /= kernel.sum()

    return kernel


def VCRfinder(
    adata,
    binsize,
    max_region,
    n_kernels=20,
    penalties=[1],
    region=None,
):
    """
    Detects variable chromatin regions (VCRs) from a anndata object containing genomic signal data
    in equally sized bins (see :ref:`scCountReads`) .

    First, a bin-to-bin correlation matrix is computed for each chromosome.

    Then, the correlation matrix is turned into a score map by convolving a number of square
    Gaussian kernels along its main diagonal. Each kernel has a sigma calculated using. Each kernel produces a 1-D score
    for each bin, which are stacked into a matrix where each row corresponds to a kernel scale and each column to a bin.

    Finally, the PELT change-point detection algorithm is applied to the score map to identify
    regions with distinct correlation patterns. This step depends on a penalty parameter that
    controls the number of detected regions.

    The function returns a pandas DataFrame containing the detected variable chromatin regions at each
    penalty. The DataFrame has columns: 'penalty', 'chrom', 'start', 'end'.

    Parameters
    ----------
    adata : anndata.AnnData
        Input anndata object with binned chromatin data. `adata.var` must contain 'chrom', 'start',
        and 'end' columns.
    binsize : int
        Size of the bins in base pairs.
    max_region : int
        Size of the largest kernel in base pairs.
    n_kernels : int, optional
        Number of Gaussian kernels to use for convolution. Default is 20.
    penalties : list of float, optional
        List of penalty values for the change-point detection algorithm. Default is [1].
    region : str, optional
        Genomic region to limit the analysis to (e.g., 'chr1:100000:200000'). Default is None.

    Returns
    -------
    output : pd.DataFrame
        Output DataFrame with detected variable chromatin regions at each penalty.
    """

    if region is not None:
        chrom, start, end = region.split(":")
        start, end = map(int, (start, end))
        adata = adata[:, (adata.var["chrom"] == chrom) & (adata.var["start"] >= start) & (adata.var["end"] <= end)]

    chroms = adata.var["chrom"].unique()
    pen_bed_df = pd.DataFrame(columns=["penalty", "chrom", "start", "end"])

    for chrom in chroms:
        adata_chrom = adata[:, adata.var["chrom"] == chrom]
        start = adata_chrom.var["start"].min()

        assert all(
            (adata_chrom.var["end"] - adata_chrom.var["start"] == binsize)
            | (adata_chrom.var["end"] - adata_chrom.var["start"] == (binsize - 1))
        ), f"Variable bin sizes detected in chromosome {chrom}"

        # Calculate sigmas
        k_factor = (max_region / binsize) ** (1 / n_kernels)
        sigmas = 0.25 * k_factor ** np.arange(1, n_kernels + 1)

        # Calculate the number of diagonals needed
        k = max(round(4.0 * np.max(sigmas)), adata_chrom.shape[1])

        ctime = time.time()
        corr = sparse_band_corr(adata_chrom.X, k=k)
        ctime = time.time() - ctime
        print(f"Band correlation with {k} diagonals calculated in {ctime:.2f} seconds")

        scores = np.zeros((len(sigmas), corr.shape[0]))

        for i, sigma in enumerate(tqdm(sigmas)):
            kernel = distance_kernel(sigma=sigma)
            k_width = kernel.shape[0]
            k_radius = k_width // 2

            # Apply kernel to diagonal
            pad_state = np.pad(corr, k_radius, mode="constant", constant_values=0.0)

            for pos in range(corr.shape[0]):
                scores[i, pos] = np.sum(pad_state[pos : pos + k_width, pos : pos + k_width] * kernel)

        # Change-point detection
        cpd = rpt.KernelCPD(kernel="rbf")

        for pen in penalties:
            try:
                bkps = cpd.fit_predict(corr, n_bkps=None, pen=pen)
            except BadSegmentationParameters:
                print(f" - No breakpoints detected for penalty {pen} in chrom {chrom}.")
                continue

            sys.stdout.write(f"Detected {len(bkps)} VCRs in chromosome {chrom} with penalty {pen}\n")

            prevs = [0, *bkps[:-1]]
            bkp_df = pd.DataFrame(
                {
                    "penalty": pen,
                    "chrom": chrom,
                    "start": [int(start + prev * binsize) for prev in prevs],
                    "end": [int(start + bkp * binsize) for bkp in bkps],
                }
            )

            pen_bed_df = pd.concat([pen_bed_df, bkp_df], ignore_index=True)

    return pen_bed_df
