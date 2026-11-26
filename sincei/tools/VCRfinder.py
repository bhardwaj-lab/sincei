from __future__ import annotations

import logging
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from typing import TYPE_CHECKING, cast

import numpy as np
import pandas as pd
import ruptures as rpt
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

log = logging.getLogger(__name__)


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


def _signal_blocks(X: CountMatrix) -> list[tuple[int, np.ndarray | sparse.csc_array]]:
    """
    Split X into blocks of contiguous features with signal (nonzero column sum).
    A feature without signal always splits two blocks.

    Parameters
    ----------
    X : scipy.sparse array or np.ndarray
        Input data matrix of shape (n_samples, n_features).

    Returns
    -------
    blocks : list of (int, np.ndarray or scipy.sparse.csc_array)
        One (start, submatrix) pair per block, in feature order. start is the
        column index of the first feature of the block in X. All features of a
        block have signal. Sparse blocks are CSC.
    """
    if not isinstance(X, np.ndarray):
        X = sparse.csc_array(X)

    active = np.flatnonzero(np.asarray(X.sum(axis=0)).ravel() != 0)
    if active.size == 0:
        return []

    breaks = np.flatnonzero(np.diff(active) > 1)
    starts = active[np.r_[0, breaks + 1]]
    ends = active[np.r_[breaks, active.size - 1]] + 1

    return [(int(s), X[:, s:e]) for s, e in zip(starts, ends, strict=True)]


def _segment_kcpd(
    X: np.ndarray | sparse.csc_array, pen: float
) -> tuple[list[int], list[float]]:
    """
    Segmentation of a block with ruptures KernelCPD (PELT) and the cosine kernel.

    The signal is the binarized block, transposed: one row per bin, one column per
    cell. The kernel of bins s and t is the cosine similarity W[s, t] of their cells
    with signal. The cost of a segment of m bins is
    C = sum_t W[t, t] - (1 / m) sum_{s, t} W[s, t], over the bins of the segment.
    The segmentation is the minimum of the sum of segment costs plus `pen` for each
    segment.

    Parameters
    ----------
    X : np.ndarray or scipy.sparse.csc_array, shape (n_samples, w)
        Signal block. All bins must have signal: the cosine kernel of a bin
        without signal is NaN.
    pen : float
        Penalty for each segment, in units of bins.

    Returns
    -------
    ends : list of int
        The end of each segment (exclusive, relative to the block start), in
        order. The last end is w.
    costs : list of float
        The cost C of each segment.
    """
    B = X != 0
    sig = np.ascontiguousarray(
        (B.toarray() if sparse.issparse(B) else B).T, dtype=np.float64
    )
    w = sig.shape[0]
    if w == 1:  # ruptures does not accept a signal of one sample
        return [1], [0.0]

    ends = rpt.KernelCPD(kernel="cosine", min_size=1, jump=1).fit(sig).predict(pen=pen)

    # All W[t, t] = 1 and sum_{s, t} W[s, t] = ||sum_t u_t||^2, with u_t the
    # normalized rows, so C = m - ||sum_t u_t||^2 / m. This does not need the w x w W.
    u = sig / np.linalg.norm(sig, axis=1, keepdims=True)
    starts = np.array([0, *ends[:-1]])
    m = np.diff([0, *ends])
    sums = np.add.reduceat(u, starts, axis=0)
    costs = m - np.einsum("ij,ij->i", sums, sums) / m
    costs = np.where(costs > 0, costs, 0.0)  # C >= 0; remove rounding errors
    return [int(e) for e in ends], [float(c) for c in costs]


def _var(adata: ad.AnnData) -> pd.DataFrame:
    return cast("pd.DataFrame", adata.var)


def _block_regions(
    block: tuple[int, np.ndarray | sparse.csc_array],
    penalty: float,
    chrom: str,
    start: int,
    binsize: int,
) -> pd.DataFrame:
    """BED rows for the segments of one signal block, which starts at bin `offset`.
    The score is the segment cost."""
    offset, X_block = block
    ends, costs = _segment_kcpd(X_block, penalty)
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
    2. Prefilter: keep the 1 counts in bins with counts in at least X_pct percent of
       the cells and in bins that are at most k bins away from a bin that passes.
    3. Split the filtered matrix into signal blocks: runs of adjacent bins with
       signal. A bin without signal always splits two blocks.
    4. Segment each block with ruptures KernelCPD (PELT) and the cosine kernel. The
       kernel of bins s and t is W[s, t], the cosine similarity of their cells with
       signal. We minimize the sum of the segment costs plus a penalty for each
       segment. The cost of a segment of m bins is
       sum_t W[t, t] - (1 / m) sum_{s, t} W[s, t] over the bins in the segment. A low
       cost means a high co-occurrence inside a segment. The penalty, in units of
       bins, controls that a block does not always separate into individual bins.
       The memory for a block of w bins is O(w * n_cells).

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
    X_pct : float, optional
        Minimum percentage of cells with signal to keep a bin. Bins below it are set
        to 0 before the blocks are made. Default is 1.0 (1%).
    k : int, optional
        Number of neighbor bins to keep to the left and to the right of bins passing
        the X_pct filter. Default is 0.
    penalty : float, optional
        Penalty for each segment, in units of bins. Higher values give fewer,
        longer regions. Default is 1.0 (the value that scFindVCRs uses).
    region : str, optional
        Genomic region to limit the analysis to (e.g., 'chr1:100000:200000').
        Default is None.
    verbose : bool, optional
        Print progress messages and warnings. Default is False.
    n_jobs : int, optional
        Number of processes to use for parallel processing. Default is
        DEFAULT_PROCESSORS (half of the available processors, at most 16).

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

    # One process pool for all chromosomes
    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        for chrom in chroms:
            sys.stdout.write(f"Processing chromosome {chrom}...\n")
            adata_chrom = adata[:, _var(adata)["chrom"] == chrom]

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
                last_size = (
                    chrom_var.loc[last_idx, "end"] - chrom_var.loc[last_idx, "start"]
                )
                if last_size >= binsize / 2:
                    _var(adata_chrom).loc[last_idx, "end"] = (
                        _var(adata_chrom).loc[last_idx, "start"] + binsize
                    )
                else:
                    # eject last row
                    sys.stderr.write(
                        f"Feature {last_idx} removed due to difference in binsize\n"
                    )
                    adata_chrom = adata_chrom[:, :-1]

            # A single bin cannot be segmented and is not a region
            if adata_chrom.shape[1] < 2:
                sys.stderr.write(
                    f"Skipping chromosome {chrom}: fewer than two bins of size "
                    f"{binsize} bp.\n"
                )
                continue

            ctime = time.time()
            X = _rare_features(cast("CountMatrix", adata_chrom.X), X_pct=X_pct, k=k)
            blocks = _signal_blocks(X)
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

            chrom_dfs = list(
                tqdm(
                    executor.map(
                        block_regions,
                        blocks,
                        chunksize=max(1, len(blocks) // (4 * n_jobs)),
                    ),
                    total=len(blocks),
                    desc=f"{chrom}: Segmenting blocks",
                    disable=not verbose,
                )
            )

            if chrom_dfs:
                chrom_df = pd.concat(chrom_dfs, ignore_index=True)
                sys.stdout.write(
                    f"Detected {len(chrom_df)} VCRs in chromosome {chrom}\n"
                )
                bed_df = pd.concat([bed_df, chrom_df], ignore_index=True)

            sys.stdout.write("\n")

    if bed_df.empty:
        log.warning(
            "No VCRs found%s. Check that the chromosomes have at least two bins, and "
            f"that bins have signal in at least {X_pct}%% of the cells (--minCellPct).",
            f" in region {region}" if region is not None else "",
        )

    return bed_df
