#!/usr/bin/env python

from __future__ import annotations

import logging
import re
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np

if TYPE_CHECKING:
    import anndata as ad
    import pandas as pd


def parse_region(region_str: str) -> tuple[str, int | None, int | None]:
    """Parse region in either CHROM or CHROM:START-END format.

    Returns a tuple (chrom, start, end). For chromosome-only input `start` and
    `end` will be returned as 0 and np.inf so every feature in the chromosome is
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
    )


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

    # If region is CHROM-only (start/end are None) select all features on that
    # chromosome and infer bounds from the selected features. Otherwise use the
    # interval mask.
    if region_start is None or region_end is None:
        region_mask = adata.var["chrom"].astype(str) == str(chrom)
    else:
        region_mask = get_region_mask(adata.var, chrom, region_start, region_end)  # ty: ignore[invalid-argument-type]

    if region_mask.sum() == 0:
        msg = f"No features found in region {region}"
        raise ValueError(msg)

    var_region: pd.DataFrame = adata.var.loc[region_mask].sort_values(["start", "end"])  # ty: ignore[unresolved-attribute]

    # If CHROM-only, infer start/end from features on that chromosome.
    if region_start is None or region_end is None:
        region_start = int(var_region["start"].min())
        region_end = int(var_region["end"].max())
    feature_idx = var_region.index
    adata_region = adata[:, feature_idx]

    X = adata_region.X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)  # ty: ignore[call-non-callable]

    if X.ndim != 2:
        msg = "Expected adata.X to be a 2D matrix."
        raise ValueError(msg)

    agg = X.sum(axis=0) if mode == "sum" else X.mean(axis=0)

    x_positions = np.arange(X.shape[1])

    # Reorder rows so higher-sum rows are on top
    row_order = np.argsort(X.sum(axis=1))[::-1]
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

    if agg.min() == 0 and agg.max() == 0:
        axes[0].set_ylim(0, 1)

    # Only min/max ticks on top Y axis
    axes[0].set_ylim(agg.min(), agg.max())
    ymin, ymax = float(np.min(agg)), float(np.max(agg))
    axes[0].set_yticks([ymin, ymax])

    # Set Y limits based on user input or data min/max or arguments.
    signal_min = signal_min or float(agg.min())
    signal_max = signal_max or float(agg.max())

    # If all values are zero, set color limits to [0, 1].
    if X_plot.min == 0 and X_plot.max == 0:
        logging.info(
            "All values in track plot are zero. Setting color limits to [0, 1]."
        )
        map_min, map_max = 0, 1

    axes[0].set_ylim(signal_min, signal_max)

    # Set heatmap color limits based on user input or data min/max or arguments.
    map_min = map_min or float(X_plot.min())
    map_max = map_max or float(X_plot.max())

    # If all values are zero, set color limits to [0, 1].
    if X_plot.min == 0 and X_plot.max == 0:
        logging.info("All values in heatmap are zero. Setting color limits to [0, 1].")
        map_min, map_max = 0, 1

    # Heatmap
    im = axes[1].imshow(
        X_plot,
        aspect="auto",
        cmap=colormap,
        interpolation="nearest",
        vmin=map_min,
        vmax=map_max,
    )

    # Keep x-axis aligned with top plot
    axes[0].set_xlim(0, X.shape[1] - 1)

    n_ticks = min(9, int(figsize[1] / 2))
    n_ticks = min(n_ticks, X.shape[1])  # Avoid more ticks than features
    tick_positions = np.linspace(0, X.shape[1] - 1, n_ticks)
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
    cax = fig.add_axes([hm_pos.x0 - 0.03, hm_pos.y0, 0.015, hm_pos.height])  # ty: ignore[no-matching-overload]
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(None)  # ty: ignore[invalid-argument-type]
    cbar.ax.yaxis.set_ticks_position("left")
    cbar.ax.yaxis.set_label_position("left")

    plt.margins(x=0.1, y=0.1)
