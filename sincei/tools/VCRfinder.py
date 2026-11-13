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


def _cooccurrence(X: CountMatrix) -> np.ndarray | sparse.csc_array:
    """
    Compute the bin-to-bin co-occurrence of a matrix X. The matrix is binarized,
    then multiplied by its transpose.

    Parameters
    ----------
    X : scipy.sparse array or np.ndarray
        Input data matrix of shape (n_samples, n_features).

    Returns
    -------
    cooc : np.ndarray or scipy.sparse.csc_array, shape (n_features, n_features)
        Symmetric co-occurrence matrix. The diagonal holds the number of cells
        with a nonzero count in each bin. Sparse when X is sparse.
    """
    if not isinstance(X, np.ndarray):
        X = sparse.csc_array(X)
    B = (X != 0).astype(np.int32)
    return B.T @ B


def _normalized_cooccurrence(X: CountMatrix) -> np.ndarray:
    """
    Compute the normalized co-occurrence W = D^-1/2 O D^-1/2 of a block, with
    O = B^T B the co-occurrence of the binarized block B and D = diag(O).

    W[s, t] is the cosine similarity of the binary columns of bins s and t, and its
    diagonal is 1. A bin without signal gets D[t, t] = 1, so its row, column and
    diagonal of W are 0.

    Parameters
    ----------
    X : scipy.sparse array or np.ndarray
        Input data matrix of shape (n_samples, n_features).

    Returns
    -------
    W : np.ndarray, shape (n_features, n_features)
        Dense, symmetric normalized co-occurrence matrix.
    """
    cooc = _cooccurrence(X)
    cooc = (cooc.toarray() if sparse.issparse(cooc) else np.asarray(cooc)).astype(
        np.float64
    )
    f = np.diag(cooc)
    D_inv_sqrt = 1.0 / np.sqrt(np.where(f > 0, f, 1.0))
    return cooc * np.outer(D_inv_sqrt, D_inv_sqrt)


def _segment_dp(W: np.ndarray, pen: float) -> list[int]:
    """
    Exact segmentation of a block by dynamic programming.

    The cost of a segment [a, b) with m = b - a bins is
    C(a, b) = sum_t W[t, t] - (1 / m) sum_{s, t} W[s, t], over the bins of the
    segment, calculated in O(1) from the 2D prefix sum of W. The segmentation is
    the exact minimum of the sum of the segment costs plus `pen` for each segment.

    Parameters
    ----------
    W : np.ndarray, shape (w, w)
        Normalized co-occurrence of the block (see `_normalized_cooccurrence`).
    pen : float
        Penalty for each segment, in units of bins.

    Returns
    -------
    ends : list of int
        The end of each segment (exclusive, relative to the block start), in
        order. The last end is w.
    """
    w = W.shape[0]
    P = np.zeros((w + 1, w + 1))
    P[1:, 1:] = W.cumsum(axis=0).cumsum(axis=1)
    q = np.r_[0.0, np.cumsum(np.diag(W))]
    F = np.zeros(w + 1)
    best = np.zeros(w + 1, dtype=np.int64)
    for b in range(1, w + 1):
        a = np.arange(b)
        inside = P[b, b] - P[a, b] - P[b, a] + P[a, a]
        total = F[:b] + (q[b] - q[a]) - inside / (b - a) + pen
        best[b] = int(np.argmin(total))
        F[b] = total[best[b]]
    ends = []
    b = w
    while b > 0:
        ends.append(b)
        b = int(best[b])
    return ends[::-1]


def _var(adata: ad.AnnData) -> pd.DataFrame:
    return cast("pd.DataFrame", adata.var)


def _segment_costs(W: np.ndarray, ends: list[int]) -> list[float]:
    """Cost of each segment of a block, sum_t W[t, t] - (1 / m) sum_{s, t} W[s, t]."""
    costs = []
    for a, b in zip([0, *ends[:-1]], ends, strict=True):
        inside = W[a:b, a:b]
        costs.append(float(np.trace(inside) - inside.sum() / (b - a)))
    return costs


def _regions_df(
    bkps: list[int],
    costs: list[float],
    chrom: str,
    start: int,
    offset: int,
    binsize: int,
) -> pd.DataFrame:
    """BED rows for the segments of one block, which starts at bin `offset`. The
    score is the segment cost."""
    prevs = [0, *bkps[:-1]]
    return pd.DataFrame({
        "chrom": chrom,
        "start": [int(start + (offset + prev) * binsize) for prev in prevs],
        "end": [int(start + (offset + bkp) * binsize) for bkp in bkps],
        "name": [f"{chrom}_VCR_{offset + bkp}" for bkp in bkps],
        "score": costs,
        "strand": "*",
    })


def _block_regions(
    block: tuple[int, np.ndarray | sparse.csc_array],
    penalty: float,
    chrom: str,
    start: int,
    binsize: int,
) -> pd.DataFrame:
    """Segment one signal block on its normalized co-occurrence."""
    offset, X_block = block

    if X_block.shape[1] == 1:
        return _regions_df([1], [0.0], chrom, start, offset, binsize)

    W = _normalized_cooccurrence(X_block)
    ends = _segment_dp(W, penalty)
    return _regions_df(ends, _segment_costs(W, ends), chrom, start, offset, binsize)


def VCRfinder(
    adata: ad.AnnData,
    binsize: int,
    penalty: float = 1.0,
    X_pct: float = 1.0,
    k: int = 0,
    region: str | None = None,
    verbose: bool = False,
    n_jobs: int = DEFAULT_PROCESSORS,
) -> pd.DataFrame:
    """
    Detects variable chromatin regions (VCRs) from an anndata object containing
    genomic signal data in equally sized bins (see :ref:`scCountReads`).

    Each chromosome is processed separately:

    1. The counts are binarized: a cell has signal in a bin if its count is not 0.
    2. Prefilter: a bin passes if at least X_pct percent of the cells have signal in
       it. The bins that pass, and the bins at most k bins from them, are kept; all
       other bins are set to 0.
    3. Signal blocks: runs of bins with signal, with at most k bins without signal
       between two bins with signal. Bins outside all blocks get no region.
    4. For each block, the co-occurrence O = B^T B of the binarized block B, and the
       normalized co-occurrence W = D^-1/2 O D^-1/2, with D = diag(O). W[s, t] is
       the cosine similarity of the cells with signal in bins s and t.
    5. Each block is segmented exactly by dynamic programming. The cost of a
       segment of m bins is sum_t W[t, t] - (1 / m) sum_{s, t} W[s, t] over its
       bins, and the segmentation minimizes the sum of the segment costs plus the
       penalty for each segment. A low cost means a high co-occurrence inside the
       segments; the penalty, in units of bins, controls the number of regions.

    The function returns a pandas DataFrame containing the detected variable
    chromatin regions. The DataFrame has columns: chrom, start, end, name, score,
    strand. The score is the segment cost (0 for a block of one bin).

    Parameters
    ----------
    adata : anndata.AnnData
        Input anndata object with binned chromatin data. `adata.var` must contain
        'chrom', 'start', and 'end' columns.
    binsize : int
        Size of the bins in base pairs.
    penalty : float, optional
        Penalty for each segment, in units of bins. Higher values give fewer,
        longer regions. Default is 1.0.
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

        # assert all(
        #     (adata_chrom.var["end"] - adata_chrom.var["start"] == binsize)
        #     | (adata_chrom.var["end"] - adata_chrom.var["start"] == (binsize - 1))
        # ), f"Variable bin sizes detected in chromosome {chrom}"

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
