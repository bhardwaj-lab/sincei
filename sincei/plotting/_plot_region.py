#!/usr/bin/env python

from __future__ import annotations

import logging
import re
from typing import TYPE_CHECKING, cast

import matplotlib.pyplot as plt
import numpy as np

if TYPE_CHECKING:
    import anndata as ad
    import pandas as pd
    from scipy.sparse import csc_matrix, csr_matrix

    # The matrix layouts `anndata.AnnData.X` holds for a count matrix.
    CountMatrix = np.ndarray | csr_matrix | csc_matrix


def parse_region(region_str: str) -> tuple[str, int | None, int | None]:
    """Parse region in either CHROM or CHROM:START-END format.

    Returns a tuple (chrom, start, end). For chromosome-only input `start` and
    `end` will be returned as None so every feature in the chromosome is
    included.
    """
    s = region_str.strip()

    # CHROM only
    if re.fullmatch(r"[^:]+", s):
        return s, None, None

    # CHROM:START-END
    match = re.fullmatch(r"([^:]+):(\d+)-(\d+)", s)
    if not match:
        msg = f"Invalid region '{region_str}'. Expected 'CHROM' or 'CHROM:START-END'"
        raise ValueError(msg)

    chrom, start_str, end_str = match.groups()
    start = int(start_str)
    end = int(end_str)
    if start > end:
        msg = f"Invalid region '{region_str}'. START must be <= END."
        raise ValueError(msg)

    return chrom, start, end


def get_region_mask(
    var_df: pd.DataFrame, chrom: str, start: int, end: int
) -> np.ndarray:
    required_cols = {"chrom", "start", "end"}
    missing = required_cols.difference(var_df.columns)
    if missing:
        msg = f"anndata.var is missing required columns: {sorted(missing)}"
        raise ValueError(msg)

    # Keep features overlapping the query interval.
    return (
        (var_df["chrom"].astype(str) == str(chrom))
        & (var_df["end"].astype(int) >= start)
        & (var_df["start"].astype(int) <= end)
    ).to_numpy()


def get_dense_columns(matrix: CountMatrix | None, col_idx: np.ndarray) -> np.ndarray:
    """Return the `col_idx` columns of a count matrix as a dense 2D array."""
    if matrix is None:
        msg = "anndata.X has no count matrix."
        raise ValueError(msg)

    selection = matrix[:, col_idx]
    # Sparse layouts carry `toarray`, dense ones do not.
    to_dense = getattr(selection, "toarray", None)
    dense = np.asarray(selection if to_dense is None else to_dense())
    if dense.ndim != 2:
        msg = "Expected adata.X to be a 2D matrix."
        raise ValueError(msg)

    return dense


# Helper function to format the unit of coordinate tick labels.
def format_genomic_tick_labels(
    min_val: int, max_val: int, n_ticks: int = 9
) -> tuple[np.ndarray, list[str], str]:
    span = max_val - min_val

    if span > 1_000_000:
        scale, unit, decimals = 1_000_000, "Mb", 1
    elif span > 5_000:
        scale, unit, decimals = 1_000, "Kb", 1
    else:
        scale, unit, decimals = 1, "bp", 0

    tick_vals = np.linspace(min_val, max_val, n_ticks)
    tick_vals = np.round(tick_vals)
    if scale == 1:
        tick_labels = [f"{int(v):,} {unit}" for v in tick_vals]
    else:
        tick_labels = [f"{v / scale:.{decimals}f} {unit}" for v in tick_vals]

    return tick_vals, tick_labels, unit


def plot_region(
    adata: ad.AnnData,
    region: str,
    mode: str = "sum",
    color: str = "red",
    colormap: str = "Reds",
    signal_min: float | None = None,
    signal_max: float | None = None,
    map_min: float | None = None,
    map_max: float | None = None,
    figsize: tuple[int, int] = (10, 6),
    dpi: int = 100,
) -> None:

    chrom, region_start, region_end = parse_region(region)
    var: pd.DataFrame = cast("pd.DataFrame", adata.var)

    # If region is CHROM-only (start/end are None) select all features on that
    # chromosome and infer bounds from the selected features. Otherwise use the
    # interval mask.
    if region_start is None or region_end is None:
        region_mask = (var["chrom"].astype(str) == str(chrom)).to_numpy()
    else:
        region_mask = get_region_mask(var, chrom, region_start, region_end)

    if region_mask.sum() == 0:
        msg = f"No features found in region {region}"
        raise ValueError(msg)

    selected = np.flatnonzero(region_mask)
    var_selected = var.iloc[selected]
    order = np.lexsort((
        var_selected["end"].to_numpy(dtype=int),
        var_selected["start"].to_numpy(dtype=int),
    ))
    col_idx = selected[order]
    var_region: pd.DataFrame = var_selected.iloc[order]

    # If CHROM-only, infer start/end from features on that chromosome.
    if region_start is None or region_end is None:
        region_start = int(var_region["start"].min())
        region_end = int(var_region["end"].max())

    X = get_dense_columns(cast("CountMatrix | None", adata.X), col_idx)

    agg = (X.sum(axis=0) if mode == "sum" else X.mean(axis=0)).ravel()
    row_totals = X.sum(axis=1).ravel()

    n_features = len(col_idx)
    x_positions = np.arange(n_features)

    # Reorder the rows so the busiest cells sit on top.
    row_order = np.ascontiguousarray(np.argsort(row_totals)[::-1])
    X_plot = X[row_order, :]

    fig, axes = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=figsize,
        sharex=True,
        gridspec_kw={"height_ratios": [1, 3]},
    )

    # Top line plot + fill
    axes[0].plot(x_positions, agg, color=color, linewidth=1.5)
    axes[0].fill_between(x_positions, agg, 0, color=color, alpha=0.9, linewidth=0)
    axes[0].set_ylabel(None)
    axes[0].tick_params(axis="x", which="both", bottom=False, labelbottom=False)

    # Y range of the track: the data range, or [0, 1] when the track is flat at
    # zero.
    ymin, ymax = float(np.min(agg)), float(np.max(agg))
    if ymin == 0.0 and ymax == 0.0:
        logging.info(
            "All values in the track plot are zero. Using a Y range of [0, 1]."
        )
        ymin, ymax = 0.0, 1.0
    axes[0].set_yticks([ymin, ymax])
    axes[0].set_ylim(
        ymin if signal_min is None else signal_min,
        ymax if signal_max is None else signal_max,
    )

    # Colour range of the heatmap, by the same rule.
    vmin, vmax = float(X_plot.min()), float(X_plot.max())
    if vmin == 0.0 and vmax == 0.0:
        logging.info(
            "All values in the heatmap are zero. Using a color range of [0, 1]."
        )
        vmin, vmax = 0.0, 1.0
    if map_min is not None:
        vmin = map_min
    if map_max is not None:
        vmax = map_max

    # Heatmap
    im = axes[1].imshow(
        X_plot,
        aspect="auto",
        cmap=colormap,
        interpolation="nearest",
        vmin=vmin,
        vmax=vmax,
    )

    # Keep x-axis aligned with top plot
    axes[0].set_xlim(0, n_features - 1)

    n_ticks = min(9, int(figsize[1] / 2))
    n_ticks = min(n_ticks, n_features)  # Avoid more ticks than features
    tick_positions = np.linspace(0, n_features - 1, n_ticks)
    _, tick_labels, _ = format_genomic_tick_labels(
        region_start, region_end, n_ticks=n_ticks
    )

    axes[1].set_xticks(tick_positions)
    axes[1].set_xticklabels(tick_labels)
    axes[1].set_xlabel(f"{chrom}:{region_start:,}-{region_end:,}", fontsize=16)
    axes[1].set_ylabel(None)
    axes[1].set_yticks([])

    # Hide all spines first
    for ax in axes:
        for spine in ax.spines.values():
            spine.set_visible(False)

    # Show and offset top-plot Y spine
    axes[0].spines["left"].set_visible(True)
    axes[0].spines["left"].set_position(("outward", 10))

    # Show and offset heatmap X spine
    axes[1].spines["bottom"].set_visible(True)
    axes[1].spines["bottom"].set_position(("outward", 10))

    # Layout first; leave room on the left for colorbar
    fig.tight_layout(rect=(0.08, 0, 1, 1))

    # Add colorbar to the LEFT of the heatmap
    hm_pos = axes[1].get_position()
    cax = fig.add_axes((hm_pos.x0 - 0.03, hm_pos.y0, 0.015, hm_pos.height))
    cbar = fig.colorbar(im, cax=cax)
    cbar.ax.yaxis.set_ticks_position("left")
    cbar.ax.yaxis.set_label_position("left")

    plt.margins(x=0.1, y=0.1)
