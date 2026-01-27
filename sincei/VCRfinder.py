import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import anndata as ad
import ruptures as rpt
from ruptures.exceptions import BadSegmentationParameters
from scipy import sparse
from tqdm import tqdm

# logs
import warnings
import logging

logger = logging.getLogger()
warnings.simplefilter(action="ignore", category=FutureWarning)


### Helper functions ###
def sparse_band_corr(X, k, chrom=None, verbose=True):
    """
    Compute only the first k diagonals of the correlation matrix of X,
    stored in banded format. Works directly on sparse matrices.

    Parameters
    ----------
    X : scipy.sparse matrix or np.ndarray
        Input data matrix of shape (n_samples, n_features).
    k : int
        Number of diagonals to compute (bandwidth).
    chrom : str or None, optional
        Chromosome name for progress bar description.
    verbose : bool, optional
        Whether to display progress bars.

    Returns
    -------
    band_corr : np.ndarray, shape (2*k+1, n_features)
        Banded correlation matrix.
    """
    n, p = X.shape
    band_corr = np.zeros((2 * k + 1, p))

    # Main diagonal is always 1
    band_corr[k, :] = 1.0

    if sparse.issparse(X):
        # Work directly with sparse matrix (convert to CSC for efficient column access)
        X_csc = X.tocsc() if not sparse.isspmatrix_csc(X) else X

        # Compute column means and norms without densifying
        col_means = np.asarray(X_csc.mean(axis=0)).ravel()

        # Compute column norms: ||x - mean||
        # = sqrt(sum(x^2) - 2*mean*sum(x) + n*mean^2)
        # = sqrt(sum(x^2) - n*mean^2)
        col_sq_sums = np.asarray(X_csc.multiply(X_csc).sum(axis=0)).ravel()
        col_norms = np.sqrt(col_sq_sums - n * col_means**2)
        col_norms[col_norms == 0] = np.inf

        # Compute correlations for each diagonal offset
        for d in tqdm(range(1, k + 1), desc=f"{chrom}: Computing banded covariance", disable=not verbose):
            n_pairs = p - d
            corr_vals = np.zeros(n_pairs)

            # Batch compute raw dot products: left_cols.T @ right_cols
            left_cols = X_csc[:, :n_pairs]
            right_cols = X_csc[:, d : d + n_pairs]

            # Compute element-wise products and sum per column pair
            # For sparse matrices, this is efficient with multiply + sum
            raw_dots = np.asarray(left_cols.multiply(right_cols).sum(axis=0)).ravel()

            # Centered dot products
            centered_dots = raw_dots - n * col_means[:n_pairs] * col_means[d : d + n_pairs]

            # Normalize to get correlations
            corr_vals = centered_dots / (col_norms[:n_pairs] * col_norms[d : d + n_pairs])

            # Store in upper and lower bands (symmetric)
            band_corr[k - d, d:] = corr_vals
            band_corr[k + d, :n_pairs] = corr_vals
    else:
        # Dense path: use efficient einsum
        # Center and normalize columns
        Xc = X - X.mean(axis=0, keepdims=True)
        norms = np.linalg.norm(Xc, axis=0)
        norms[norms == 0] = np.inf
        Xn = Xc / norms

        # Off-diagonals
        for d in tqdm(range(1, k + 1), desc=f"{chrom}: Computing banded covariance", disable=not verbose):
            corr_vals = np.einsum(
                "ij,ij->j",
                Xn[:, : p - d],
                Xn[:, d:],
                optimize=True,
            )
            band_corr[k - d, d:] = corr_vals
            band_corr[k + d, : p - d] = corr_vals

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
    verbose=False,
    n_threads=None,
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
    verbose : bool, optional
        Print additional messages.

    Returns
    -------
    output : pd.DataFrame
        Output DataFrame with detected variable chromatin regions at each penalty.
    """
    adata.var[["start", "end"]] = adata.var[["start", "end"]].apply(pd.to_numeric)

    if n_threads is None or n_threads < 1:
        n_threads = 1

    if region is not None:
        region = region.split(":")
        if len(region) == 3:
            chrom, start, end = region[0], int(region[1]), int(region[2])
        elif len(region) == 2:
            chrom, start = region[0], int(region[1])
            end = np.max(adata.var.loc[adata.var.chrom == chrom, "end"])
        else:
            chrom = region[0]
            start = np.min(adata.var.loc[adata.var.chrom == chrom, "start"])
            end = np.max(adata.var.loc[adata.var.chrom == chrom, "end"])

        adata = adata[:, (adata.var["chrom"] == chrom) & (adata.var["start"] >= start) & (adata.var["end"] <= end)]

    chroms = adata.var["chrom"].unique()
    pen_bed_df = pd.DataFrame(columns=["chrom", "start", "end", "name", "score", "strand"])

    for chrom in chroms:
        sys.stdout.write(f"Processing chromosome {chrom}...\n")
        adata_chrom = adata[:, adata.var["chrom"] == chrom]

        if adata_chrom.shape[1] == 1:
            sys.stderr.write(f"Skipping chromosome {chrom} with only one bin.\n")
            bkp_df = pd.DataFrame(
                {
                    "chrom": [chrom] * len(penalties),
                    "start": [int(adata_chrom.var["start"].values[0])] * len(penalties),
                    "end": [int(adata_chrom.var["end"].values[0])] * len(penalties),
                    "name": [f"pen-{pen}_brkpoint-1" for pen in penalties],
                    "score": penalties,
                    "strand": ["*"] * len(penalties),
                }
            )
            pen_bed_df = pd.concat([pen_bed_df, bkp_df], ignore_index=True)
            continue

        start = adata_chrom.var["start"].min()
        # Sort the bins
        adata_chrom = adata_chrom[:, adata_chrom.var.sort_values(by=["start", "end"], axis=0).index]

        # rows where end - start != binsize
        mask = (adata_chrom.var["end"] - adata_chrom.var["start"]) != binsize

        # if the last row binsize is >= 1/2 of binsize, update its 'end' value
        last_idx = adata_chrom.var.index[-1]

        if mask.loc[last_idx]:
            if adata_chrom.var.loc[last_idx, "end"] - adata_chrom.var.loc[last_idx, "start"] >= (binsize / 2):
                adata_chrom.var.loc[last_idx, "end"] = adata_chrom.var.loc[last_idx, "start"] + binsize
            else:
                # eject last row
                sys.stderr.write(f"Feature {last_idx} removed due to difference in binsize\n")
                adata_chrom = adata_chrom[:, : adata_chrom.var.index[-2]]

        assert all(
            (adata_chrom.var["end"] - adata_chrom.var["start"] == binsize)
            | (adata_chrom.var["end"] - adata_chrom.var["start"] == (binsize - 1))
        ), f"Variable bin sizes detected in chromosome {chrom}"

        reg_to_consider = min(max_region, max(int(adata_chrom.shape[1] * binsize / n_kernels), binsize))

        # Calculate sigmas
        k_factor = (reg_to_consider / binsize) ** (1 / n_kernels)
        sigmas = 0.25 * k_factor ** np.arange(1, n_kernels + 1)

        # Calculate the number of diagonals needed (capped at number of bins)
        k = min(round(4.0 * np.max(sigmas)), adata_chrom.shape[1])

        # Get bin-bin correlations
        ctime = time.time()
        corr = sparse_band_corr(adata_chrom.X, k=k, chrom=chrom, verbose=verbose)
        p = adata_chrom.shape[1]  # number of bins
        ctime = time.time() - ctime
        if verbose:
            sys.stdout.write(
                f"Chromosome {chrom}: Band correlation with {k} diagonals calculated in {ctime:.2f} seconds\n"
            )

        scores = np.zeros((p, len(sigmas)))

        def _score_sigma(args):
            i, sigma = args
            radius = min(k, int(4.0 * sigma))
            kernel = distance_kernel(sigma=sigma, radius=radius)
            k_width = kernel.shape[0]
            k_radius = k_width // 2

            score_row = np.zeros(p)

            # Decompose 2D kernel convolution into 1D convolutions per diagonal
            # score[pos] = sum over (di, dj) of kernel[di+r, dj+r] * corr[pos+di, pos+dj]
            # Group by d = dj - di (the diagonal offset in corr)
            for d in range(-min(k_radius, k), min(k_radius, k) + 1):
                # Extract diagonal d from kernel: elements where col - row = d
                kernel_diag = np.diag(kernel, k=d)  # shape: (k_width - |d|,)

                # Get correlation values for diagonal offset d from banded storage
                # band_corr[k - d, :] contains corr[i, i+d] for d >= 0
                # band_corr[k - d, :] contains corr[i-d, i] for d < 0
                corr_diag = corr[k - d, :]  # shape: (p,)

                # Convolve: this computes sum of kernel_diag[j] * corr_diag[pos + j - offset]
                # We need to account for the offset in the kernel diagonal
                conv = np.convolve(corr_diag, kernel_diag[::-1], mode="same")
                score_row += conv

            return i, score_row

        if n_threads > 1:
            with ThreadPoolExecutor(max_workers=n_threads) as executor:
                for i, score_row in tqdm(
                    executor.map(_score_sigma, list(enumerate(sigmas))),
                    total=len(sigmas),
                    desc=f"{chrom}: Calculating score matrix",
                    disable=not verbose,
                ):
                    scores[:, i] = score_row
        else:
            for i, sigma in tqdm(
                list(enumerate(sigmas)), desc=f"{chrom}: Calculating score matrix", disable=not verbose
            ):
                _, score_row = _score_sigma((i, sigma))
                scores[:, i] = score_row

        # Use Pelt with L2 cost - O(n) memory
        # scores has shape (p, n_kernels) - each position is a feature vector
        algo = rpt.KernelCPD(kernel="linear", min_size=1, jump=1).fit(scores)

        def _predict_penalty(pen):
            """Predict breakpoints for a single penalty value."""
            try:
                bkps = algo.predict(pen=pen)
            except BadSegmentationParameters:
                sys.stderr.write(f" - No breakpoints detected for penalty {pen} in chrom {chrom}.\n")
                return None

            sys.stdout.write(f"Detected {len(bkps)} VCRs in chromosome {chrom} with penalty {pen}\n")

            prevs = [0, *bkps[:-1]]
            bkp_df = pd.DataFrame(
                {
                    "chrom": chrom,
                    "start": [int(start + prev * 2000) for prev in prevs],
                    "end": [int(start + bkp * 2000) for bkp in bkps],
                    "name": [f"pen-{pen}_brkpoint-{bkp}" for bkp in bkps],
                    "score": pen,
                    "strand": "*",
                }
            )
            return bkp_df

        if n_threads > 1:
            with ThreadPoolExecutor(max_workers=n_threads) as executor:
                for bkp_df in tqdm(
                    executor.map(_predict_penalty, penalties),
                    total=len(penalties),
                    desc=f"{chrom}: Change-point detection",
                    disable=not verbose,
                ):
                    if bkp_df is not None:
                        pen_bed_df = pd.concat([pen_bed_df, bkp_df], ignore_index=True)
        else:
            for pen in tqdm(penalties, desc=f"{chrom}: Change-point detection", disable=not verbose):
                bkp_df = _predict_penalty(pen)
                if bkp_df is not None:
                    pen_bed_df = pd.concat([pen_bed_df, bkp_df], ignore_index=True)

        sys.stdout.write("\n")

    return pen_bed_df
