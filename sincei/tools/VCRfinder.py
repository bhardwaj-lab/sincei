from __future__ import annotations

import sys
import time
import warnings
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from typing import TYPE_CHECKING, cast

import numpy as np
import pandas as pd
import ruptures as rpt
from ruptures.exceptions import BadSegmentationParameters
from scipy import sparse
from tqdm import tqdm

if TYPE_CHECKING:
    from collections.abc import Sequence

    import anndata as ad

    # The matrix layouts `anndata.AnnData.X` holds for a count matrix.
    CountMatrix = (
        np.ndarray
        | sparse.csr_matrix
        | sparse.csc_matrix
        | sparse.csr_array
        | sparse.csc_array
    )

warnings.simplefilter(action="ignore", category=FutureWarning)


### Helper functions ###
def sparse_band_corr(
    X: CountMatrix,
    k: int,
    chrom: str | None = None,
    verbose: bool = True,
) -> np.ndarray:
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
        X_csc = sparse.csc_matrix(X)

        # Compute column means and norms without densifying
        col_means = np.asarray(X_csc.mean(axis=0)).ravel()

        # Compute column norms: ||x - mean||
        # = sqrt(sum(x^2) - 2*mean*sum(x) + n*mean^2)
        # = sqrt(sum(x^2) - n*mean^2)
        col_sq_sums = np.asarray(X_csc.multiply(X_csc).sum(axis=0)).ravel()
        col_norms = np.sqrt(col_sq_sums - n * col_means**2)
        col_norms[col_norms == 0] = np.inf

        # Compute correlations for each diagonal offset
        for d in tqdm(
            range(1, k + 1),
            desc=f"{chrom}: Computing banded covariance",
            disable=not verbose,
        ):
            n_pairs = p - d

            # Batch compute raw dot products: left_cols.T @ right_cols
            left_cols = X_csc[:, :n_pairs]
            right_cols = X_csc[:, d : d + n_pairs]

            # Compute element-wise products and sum per column pair
            # For sparse matrices, this is efficient with multiply + sum
            raw_dots = np.asarray(left_cols.multiply(right_cols).sum(axis=0)).ravel()

            # Centered dot products
            centered_dots = (
                raw_dots - n * col_means[:n_pairs] * col_means[d : d + n_pairs]
            )

            # Normalize to get correlations
            corr_vals = centered_dots / (
                col_norms[:n_pairs] * col_norms[d : d + n_pairs]
            )

            # Store in upper and lower bands (symmetric)
            band_corr[k - d, d:] = corr_vals
            band_corr[k + d, :n_pairs] = corr_vals
    else:
        # Dense path: use efficient einsum
        # Center and normalize columns
        dense = np.asarray(X)
        Xc = dense - dense.mean(axis=0, keepdims=True)
        norms = np.linalg.norm(Xc, axis=0)
        norms[norms == 0] = np.inf
        Xn = Xc / norms

        # Off-diagonals
        for d in tqdm(
            range(1, k + 1),
            desc=f"{chrom}: Computing banded covariance",
            disable=not verbose,
        ):
            corr_vals = np.einsum(
                "ij,ij->j",
                Xn[:, : p - d],
                Xn[:, d:],
                optimize=True,
            )
            band_corr[k - d, d:] = corr_vals
            band_corr[k + d, : p - d] = corr_vals

    return band_corr


def distance_kernel(
    sigma: float, truncate: float = 4.0, radius: int | None = None
) -> np.ndarray:
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


def _var(adata: ad.AnnData) -> pd.DataFrame:
    return cast("pd.DataFrame", adata.var)


def _score_sigma(corr: np.ndarray, k: int, sigma: float) -> np.ndarray:
    """Score every bin of a banded correlation matrix with one Gaussian kernel."""
    radius = min(k, int(4.0 * sigma))
    kernel = distance_kernel(sigma=sigma, radius=radius)
    k_width = kernel.shape[0]
    k_radius = k_width // 2

    score_row = np.zeros(corr.shape[1])

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

        # Convolve: this computes sum of
        # kernel_diag[j] * corr_diag[pos + j - offset]
        # We need to account for the offset in the kernel diagonal
        conv = np.convolve(corr_diag, kernel_diag[::-1], mode="same")
        score_row += conv

    return score_row


def _penalty_regions(
    algo: rpt.KernelCPD, pen: float, chrom: str, start: int, binsize: int
) -> pd.DataFrame | None:
    """Predict breakpoints for a single penalty value."""
    try:
        bkps = algo.predict(pen=pen)
    except BadSegmentationParameters:
        sys.stderr.write(
            f" - No breakpoints detected for penalty {pen} in chrom {chrom}.\n"
        )
        return None

    sys.stdout.write(
        f"Detected {len(bkps)} VCRs in chromosome {chrom} with penalty {pen}\n"
    )

    prevs = [0, *bkps[:-1]]
    return pd.DataFrame({
        "chrom": chrom,
        "start": [int(start + prev * binsize) for prev in prevs],
        "end": [int(start + bkp * binsize) for bkp in bkps],
        "name": [f"{chrom}_VCR_{bkp}_pen{pen}" for bkp in bkps],
        "score": pen,
        "strand": "*",
    })


def VCRfinder(
    adata: ad.AnnData,
    binsize: int,
    max_region: int,
    n_kernels: int = 20,
    penalties: Sequence[float] = (1,),
    region: str | None = None,
    verbose: bool = False,
    n_threads: int = 1,
) -> pd.DataFrame:
    """
    Detects variable chromatin regions (VCRs) from a anndata object containing
    genomic signal data in equally sized bins (see :ref:`scCountReads`) .

    First, a bin-to-bin correlation matrix is computed for each chromosome.

    Then, the correlation matrix is turned into a score map by convolving a number
    of square Gaussian kernels along its main diagonal. Each kernel has a sigma
    calculated using a maximum region size to consider. Each kernel produces a 1-D
    score for each bin, which are stacked into a matrix where each row corresponds
    to a kernel scale and each column to a bin.

    Finally, the PELT change-point detection algorithm is applied to the score map
    to identify regions with distinct correlation patterns. This step depends on a
    penalty parameter that controls the number of detected regions.

    The function returns a pandas DataFrame containing the detected variable
    chromatin regions at each penalty. The DataFrame has columns: chrom, start, end,
    name, score, strand.

    Parameters
    ----------
    adata : anndata.AnnData
        Input anndata object with binned chromatin data. `adata.var` must contain
        'chrom', 'start', and 'end' columns.
    binsize : int
        Size of the bins in base pairs.
    max_region : int
        Size of the largest kernel in base pairs.
    n_kernels : int, optional
        Number of Gaussian kernels to use for convolution. Default is 20.
    penalties : list of float, optional
        List of penalty values for the change-point detection algorithm. Default is
        [1].
    region : str, optional
        Genomic region to limit the analysis to (e.g., 'chr1:100000:200000').
        Default is None.
    verbose : bool, optional
        Print progress messages and warnings. Default is False.
    n_threads : int, optional
        Number of threads to use for parallel processing, by default 1.

    Returns
    -------
    output : pd.DataFrame
        Output DataFrame with detected variable chromatin regions at each penalty.
    """
    _var(adata)[["start", "end"]] = _var(adata)[["start", "end"]].apply(pd.to_numeric)

    if region is not None:
        parts = region.split(":")
        region_chrom = parts[0]
        on_chrom = _var(adata).loc[_var(adata)["chrom"] == region_chrom]
        if len(parts) == 3:
            region_start, region_end = int(parts[1]), int(parts[2])
        elif len(parts) == 2:
            region_start, region_end = int(parts[1]), np.max(on_chrom["end"])
        else:
            region_start = np.min(on_chrom["start"])
            region_end = np.max(on_chrom["end"])

        adata = adata[
            :,
            (_var(adata)["chrom"] == region_chrom)
            & (_var(adata)["start"] >= region_start)
            & (_var(adata)["end"] <= region_end),
        ]

    chroms = _var(adata)["chrom"].unique()
    pen_bed_df = pd.DataFrame(
        columns=pd.Index(["chrom", "start", "end", "name", "score", "strand"])
    )

    for chrom in chroms:
        sys.stdout.write(f"Processing chromosome {chrom}...\n")
        adata_chrom = adata[:, _var(adata)["chrom"] == chrom]

        if adata_chrom.shape[1] == 1:
            sys.stderr.write(f"Skipping chromosome {chrom} with only one bin.\n")
            bkp_df = pd.DataFrame({
                "chrom": [chrom] * len(penalties),
                "start": [int(_var(adata_chrom)["start"].to_numpy()[0])]
                * len(penalties),
                "end": [int(_var(adata_chrom)["end"].to_numpy()[0])] * len(penalties),
                "name": [f"pen-{pen}_brkpoint-1" for pen in penalties],
                "score": list(penalties),
                "strand": ["*"] * len(penalties),
            })
            pen_bed_df = pd.concat([pen_bed_df, bkp_df], ignore_index=True)
            continue

        start = _var(adata_chrom)["start"].min()
        # Sort the bins
        adata_chrom = adata_chrom[
            :, _var(adata_chrom).sort_values(by=["start", "end"], axis=0).index
        ]

        # rows where end - start != binsize
        chrom_var = _var(adata_chrom)
        mask = (chrom_var["end"] - chrom_var["start"]) != binsize

        # if the last row binsize is >= 1/2 of binsize, update its 'end' value
        last_idx = chrom_var.index[-1]

        if mask.loc[last_idx]:
            if chrom_var.loc[last_idx, "end"] - chrom_var.loc[last_idx, "start"] >= (
                binsize / 2
            ):
                _var(adata_chrom).loc[last_idx, "end"] = (
                    _var(adata_chrom).loc[last_idx, "start"] + binsize
                )
            else:
                # eject last row
                sys.stderr.write(
                    f"Feature {last_idx} removed due to difference in binsize\n"
                )
                adata_chrom = adata_chrom[:, : _var(adata_chrom).index[-2]]

        # assert all(
        #     (adata_chrom.var["end"] - adata_chrom.var["start"] == binsize)
        #     | (adata_chrom.var["end"] - adata_chrom.var["start"] == (binsize - 1))
        # ), f"Variable bin sizes detected in chromosome {chrom}"

        reg_to_consider = min(
            max_region, max(int(adata_chrom.shape[1] * binsize / n_kernels), binsize)
        )

        # Calculate sigmas
        k_factor = (reg_to_consider / binsize) ** (1 / n_kernels)
        sigmas = 0.25 * k_factor ** np.arange(1, n_kernels + 1)

        # Calculate the number of diagonals needed (capped at number of bins)
        k = min(round(4.0 * np.max(sigmas)), adata_chrom.shape[1])

        # Get bin-bin correlations
        ctime = time.time()
        corr = sparse_band_corr(
            cast("CountMatrix", adata_chrom.X), k=k, chrom=chrom, verbose=verbose
        )
        p = adata_chrom.shape[1]  # number of bins
        ctime = time.time() - ctime
        if verbose:
            sys.stdout.write(
                f"Chromosome {chrom}: Banded correlation with {k} diagonals "
                f"calculated in {ctime:.2f} seconds\n"
            )

        scores = np.zeros((p, len(sigmas)))
        score_sigma = partial(_score_sigma, corr, k)

        if n_threads > 1:
            with ThreadPoolExecutor(max_workers=n_threads) as executor:
                for i, score_row in tqdm(
                    enumerate(executor.map(score_sigma, sigmas)),
                    total=len(sigmas),
                    desc=f"{chrom}: Calculating score matrix",
                    disable=not verbose,
                ):
                    scores[:, i] = score_row
        else:
            for i, sigma in tqdm(
                list(enumerate(sigmas)),
                desc=f"{chrom}: Calculating score matrix",
                disable=not verbose,
            ):
                scores[:, i] = score_sigma(sigma)

        # Use Pelt with L2 cost - O(n) memory
        # scores has shape (p, n_kernels) - each position is a feature vector
        algo = rpt.KernelCPD(kernel="linear", min_size=1, jump=1).fit(scores)
        penalty_regions = partial(
            _penalty_regions, algo, chrom=chrom, start=start, binsize=binsize
        )

        if n_threads > 1:
            with ThreadPoolExecutor(max_workers=n_threads) as executor:
                for bkp_df in tqdm(
                    executor.map(penalty_regions, penalties),
                    total=len(penalties),
                    desc=f"{chrom}: Change-point detection",
                    disable=not verbose,
                ):
                    if bkp_df is not None:
                        pen_bed_df = pd.concat([pen_bed_df, bkp_df], ignore_index=True)
        else:
            for pen in tqdm(
                penalties, desc=f"{chrom}: Change-point detection", disable=not verbose
            ):
                bkp_df = penalty_regions(pen)
                if bkp_df is not None:
                    pen_bed_df = pd.concat([pen_bed_df, bkp_df], ignore_index=True)

        sys.stdout.write("\n")

    return pen_bed_df
