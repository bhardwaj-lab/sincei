from __future__ import annotations

import sys
import time
import warnings
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from typing import TYPE_CHECKING, cast

import numpy as np
import pandas as pd
from scipy import sparse
from tqdm import tqdm

from sincei.cli._common_args import DEFAULT_PROCESSORS

if TYPE_CHECKING:
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
def _with_neighbors(keep: np.ndarray, k: int) -> np.ndarray:
    """Extend a boolean feature mask by k features to the left and to the right."""
    if k <= 0:
        return keep
    p = keep.size
    csum = np.r_[0, np.cumsum(keep)]
    idx = np.arange(p)
    return (csum[np.minimum(idx + k + 1, p)] - csum[np.maximum(idx - k, 0)]) > 0


def _rare_features(
    X: CountMatrix, X_pct: float = 1.0, k: int = 0
) -> np.ndarray | sparse.csr_array | sparse.csc_array:
    """
    Set to 0 the columns of features with nonzero counts in less than X_pct percent
    of n_samples, unless they are at most k features from a feature that passes.

    Parameters
    ----------
    X : scipy.sparse array or np.ndarray
        Input data matrix of shape (n_samples, n_features).
    X_pct : float, optional
        Minimum percentage of cells with a nonzero count to keep a feature.
        Default is 1.0 (1%).
    k : int, optional
        Number of neighbors to the left and to the right of each kept feature
        that also keep their counts. Default is 0 (no neighbors).

    Returns
    -------
    X_filtered : np.ndarray or scipy.sparse array, shape (n_samples, n_features)
        Copy of X with rare features set to 0. A CSC input is returned as
        csc_array; other sparse inputs are returned as csr_array.
    """
    n = X.shape[0]

    if not isinstance(X, np.ndarray):
        if X.format == "csc":
            Xs = sparse.csc_array(X, copy=True)
        else:
            Xs = sparse.csr_array(X, copy=True)
        Xs.eliminate_zeros()
        if Xs.format == "csc":
            nnz = np.diff(Xs.indptr)
            rare = ~_with_neighbors(nnz >= n * X_pct / 100, k)
            Xs.data[np.repeat(rare, nnz)] = 0
        else:
            nnz = np.bincount(Xs.indices, minlength=Xs.shape[1])
            rare = ~_with_neighbors(nnz >= n * X_pct / 100, k)
            Xs.data[rare[Xs.indices]] = 0
        Xs.eliminate_zeros()
        return Xs

    dense = np.array(X, copy=True)
    rare = ~_with_neighbors(np.count_nonzero(dense, axis=0) >= n * X_pct / 100, k)
    dense[:, rare] = 0
    return dense


def _signal_blocks(
    X: CountMatrix, k: int = 0
) -> list[tuple[int, np.ndarray | sparse.csc_array]]:
    """
    Split X into blocks of contiguous features with signal (nonzero column sum).

    Parameters
    ----------
    X : scipy.sparse array or np.ndarray
        Input data matrix of shape (n_samples, n_features).
    k : int, optional
        Slack: the maximum number of features with a zero column sum between two
        features with signal in the same block. These zero columns are included
        in the block. Default is 0 (only adjacent features).

    Returns
    -------
    blocks : list of (int, np.ndarray or scipy.sparse.csc_array)
        One (start, submatrix) pair per block, in feature order. start is the
        column index of the first feature of the block in X. Each block starts
        and ends with a feature with signal. Sparse blocks are CSC.
    """
    if not isinstance(X, np.ndarray):
        X = sparse.csc_array(X)

    active = np.flatnonzero(np.asarray(X.sum(axis=0)).ravel() != 0)
    if active.size == 0:
        return []

    breaks = np.flatnonzero(np.diff(active) > k + 1)
    starts = active[np.r_[0, breaks + 1]]
    ends = active[np.r_[breaks, active.size - 1]] + 1

    return [(int(s), X[:, s:e]) for s, e in zip(starts, ends, strict=True)]


def _band_cooccurrence(X: CountMatrix, L: int) -> np.ndarray:
    """
    Band of the normalized co-occurrence W = D^-1/2 O D^-1/2 of a block, with
    O = B^T B the co-occurrence of the binarized block B and D = diag(O).
    W[s, t] is the cosine similarity of the cells with signal in bins s and t. A
    bin without signal gets D[t, t] = 1, so its row and column of W are 0. The
    full W is not calculated.

    Parameters
    ----------
    X : scipy.sparse array or np.ndarray
        Input data matrix of shape (n_samples, n_features).
    L : int
        Width of the band: the maximum segment length.

    Returns
    -------
    band : np.ndarray, shape (n_features, L)
        band[s, d] = W[s, s + d], and 0 where s + d is outside the block.
    """
    B = sparse.csc_array(X != 0, dtype=np.float64)
    f = np.asarray(B.sum(axis=0)).ravel()
    B = sparse.csc_array(B @ sparse.diags_array(1.0 / np.sqrt(np.where(f > 0, f, 1.0))))
    w = B.shape[1]
    band = np.zeros((w, L))
    for c0 in range(0, w, L):
        c1 = min(c0 + L, w)
        G = (B[:, c0:c1].T @ B[:, c0 : min(c1 + L, w)]).toarray()
        for s in range(c0, c1):
            n = min(L, w - s)
            band[s, :n] = G[s - c0, s - c0 : s - c0 + n]
    return band


def _segment_dp(band: np.ndarray, pen: float) -> tuple[list[int], list[float]]:
    """
    Segmentation of a block by dynamic programming on the band of W.

    The cost of a segment [a, b) with m = b - a bins is
    C(a, b) = sum_t W[t, t] - (1 / m) sum_{s, t} W[s, t], over the bins of the
    segment. The segmentation is the minimum of the sum of segment costs plus `pen`
    for each segment, with segments of at most L bins (the width of the band). It
    is exact if L is larger than the number of bins in the block. If two segmentations
    have the same minimum, the one with the longer segment is used.

    Parameters
    ----------
    band : np.ndarray, shape (w, L)
        Band of the normalized co-occurrence (see `_band_cooccurrence`).
    pen : float
        Penalty for each segment, in units of bins.

    Returns
    -------
    ends : list of int
        The end of each segment (exclusive, relative to the block start), in
        order. The last end is w.
    costs : list of float
        The cost C(a, b) of each segment.
    """
    w, L = band.shape
    rowcum = np.cumsum(band, axis=1)
    q = np.r_[0.0, np.cumsum(band[:, 0])]
    F = np.zeros(w + 1)
    best = np.zeros(w + 1, dtype=np.int64)
    cost = np.zeros(w + 1)
    for b in range(1, w + 1):
        a = np.arange(max(b - L, 0), b)
        row_sums = 2 * rowcum[a, b - 1 - a] - band[a, 0]
        inside = np.cumsum(row_sums[::-1])[::-1]
        c = q[b] - q[a] - inside / (b - a)
        total = F[a] + c + pen
        i = int(np.argmax(total <= total.min() + 1e-9))
        best[b], cost[b], F[b] = a[i], c[i], total[i]
    ends, costs = [], []
    b = w
    while b > 0:
        ends.append(b)
        costs.append(float(cost[b]))
        b = int(best[b])
    return ends[::-1], costs[::-1]


def _var(adata: ad.AnnData) -> pd.DataFrame:
    return cast("pd.DataFrame", adata.var)


def _block_regions(
    block: tuple[int, np.ndarray | sparse.csc_array],
    penalty: float,
    max_length: int,
    chrom: str,
    start: int,
    binsize: int,
) -> pd.DataFrame:
    """BED rows for the segments of one signal block, which starts at bin `offset`.
    The score is the segment cost."""
    offset, X_block = block
    band = _band_cooccurrence(X_block, min(X_block.shape[1], max_length))
    ends, costs = _segment_dp(band, penalty)
    starts = [0, *ends[:-1]]
    return pd.DataFrame({
        "chrom": chrom,
        "start": [int(start + (offset + a) * binsize) for a in starts],
        "end": [int(start + (offset + b) * binsize) for b in ends],
        "name": [f"{chrom}_VCR_{offset + b}" for b in ends],
        "score": costs,
        "strand": "*",
    })


def VCRfinder(
    adata: ad.AnnData,
    X_pct: float = 1.0,
    k: int = 0,
    max_length: int = 2000,
    penalty: float = 1.0,
    region: str | None = None,
    verbose: bool = False,
    n_jobs: int = DEFAULT_PROCESSORS,
) -> pd.DataFrame:
    """
    Detects variable chromatin regions (VCRs) from an anndata object containing
    genomic signal data in equally sized bins (see :ref:`scCountReads`).

    Each chromosome is processed separately:

    1. Binarize the count matrix: set all nonzero counts to 1.
    2. Prefilter: keep the 1 counts in bins with at counts in at least X_pct cells and
       in bins that are at most k bins away from a bin that passes.
    3. Split the filtered matrix into signal blocks, runs of bins with signal with at
       most k bins without signal between.
    4. For each block, the compute the normalized co-occurrence the binarized block B,
       W = D^-1/2 O D^-1/2, with O the co-occurrence matrix and D = diag(O). Only the
       band of W with |s - t| < L is calculated, with L = max_length. This avoids the
       quadratic blowup of memory needed to store the whole co-occurrence matrix of a
       large block.
    5. Each block is segmented by dynamic programming. We minimize the sum of the cost
       of splitting a block into segments plus a penalty. The cost of a segment of m
       bins is sum_t W[t, t] - (1 / m) sum_{s, t} W[s, t] over the bins in the segment.
       Segments can be at most L bins long. A low cost means a high co-occurrence inside
       a segment. The penalty, in units of bins, controls that a block does not always
       separate into individual bins.

    The function returns a pandas DataFrame containing the detected variable
    chromatin regions. The DataFrame has columns: chrom, start, end, name, score,
    strand. The score is the segment cost (0 for a block of one bin).

    Parameters
    ----------
    adata : anndata.AnnData
        Input anndata object with binned chromatin data. `adata.var` must contain
        'chrom', 'start', and 'end' columns. The bin size is the most common
        end - start; more than 90% of the bins must have it, else a ValueError is
        raised.
    penalty : float, optional
        Penalty for each segment, in units of bins. Higher values give fewer,
        longer regions. Default is 1.0.
    max_length : int, optional
        Maximum number of bins in a region. Only the band of W of this width is
        calculated, so the memory for a block of w bins is O(w * max_length). The
        segmentation is exact for blocks of at most max_length bins. Default is 2000.
    X_pct : float, optional
        Minimum percentage of cells with signal to keep a bin. Bins below it are set
        to 0 before the blocks are made. Default is 1.0 (1%).
    k : int, optional
        Number of neighbor bins to the left and to the right of each bin that
        passes X_pct that are also kept in the prefilter, and the maximum number of
        bins without signal between two bins with signal in the same block.
        Default is 0.
    region : str, optional
        Genomic region to limit the analysis to (e.g., 'chr1:100000:200000').
        Default is None.
    verbose : bool, optional
        Print progress messages and warnings. Default is False.
    n_jobs : int, optional
        Number of threads to use for parallel processing. Default is DEFAULT_PROCESSORS
        (half of the available processors, at most 16).

    Returns
    -------
    output : pd.DataFrame
        Output DataFrame with the detected variable chromatin regions.
    """
    _var(adata)[["start", "end"]] = _var(adata)[["start", "end"]].apply(pd.to_numeric)

    # Detect bin size from AnnData
    sizes = (_var(adata)["end"] - _var(adata)["start"]).to_numpy()
    values, counts = np.unique(sizes, return_counts=True)
    binsize = int(values[np.argmax(counts)])
    fraction = counts.max() / len(sizes)
    if fraction <= 0.9:
        msg = (
            f"The bins in the input have varying sizes: {1 - fraction:.1%} of bins "
            f"are not of the most common size, {binsize} bp. More than 90% of the bins "
            f"must have the same size. You can make a matrix of equal-size bins with "
            f"scCountReads bins."
        )
        raise ValueError(msg)
    if verbose:
        sys.stdout.write(
            f"Bin size: {binsize} bp ({fraction:.1%} of bins are this length.)\n"
        )

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
    bed_df = pd.DataFrame(
        columns=pd.Index(["chrom", "start", "end", "name", "score", "strand"])
    )

    for chrom in chroms:
        sys.stdout.write(f"Processing chromosome {chrom}...\n")
        adata_chrom = adata[:, _var(adata)["chrom"] == chrom]

        if adata_chrom.shape[1] == 1:
            sys.stderr.write(f"Skipping chromosome {chrom} with only one bin.\n")
            bkp_df = pd.DataFrame({
                "chrom": [chrom],
                "start": [int(_var(adata_chrom)["start"].to_numpy()[0])],
                "end": [int(_var(adata_chrom)["end"].to_numpy()[0])],
                "name": [f"{chrom}_VCR_1"],
                "score": [0.0],
                "strand": ["*"],
            })
            bed_df = pd.concat([bed_df, bkp_df], ignore_index=True)
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

        ctime = time.time()
        X = _rare_features(cast("CountMatrix", adata_chrom.X), X_pct=X_pct, k=k)
        blocks = _signal_blocks(X, k=k)
        if verbose:
            sys.stdout.write(
                f"Chromosome {chrom}: {len(blocks)} signal blocks found in "
                f"{time.time() - ctime:.2f} seconds\n"
            )

        block_regions = partial(
            _block_regions,
            penalty=penalty,
            max_length=max_length,
            chrom=chrom,
            start=start,
            binsize=binsize,
        )

        with ThreadPoolExecutor(max_workers=n_jobs) as executor:
            chrom_dfs = list(
                tqdm(
                    executor.map(block_regions, blocks),
                    total=len(blocks),
                    desc=f"{chrom}: Segmenting blocks",
                    disable=not verbose,
                )
            )

        if chrom_dfs:
            chrom_df = pd.concat(chrom_dfs, ignore_index=True)
            sys.stdout.write(f"Detected {len(chrom_df)} VCRs in chromosome {chrom}\n")
            bed_df = pd.concat([bed_df, chrom_df], ignore_index=True)

        sys.stdout.write("\n")

    return bed_df
