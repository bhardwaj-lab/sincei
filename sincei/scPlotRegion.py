#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re
import sys

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np

from sincei import ParserCommon


def get_args():
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group("Plot Options")

    group.add_argument(
        "-r",
        "--region",
        required=True,
        help="Genomic region in CHROM:START-END format",
    )

    group.add_argument(
        "-m",
        "--mode",
        required=True,
        choices=["sum", "mean"],
        default="sum",
        help="Aggregation mode for the top subplot",
    )

    group.add_argument(
        "--summarize",
        action="store_true",
        help="Sum contiguous columns to reduce matrix width for heatmap plotting (default: False)",
    )

    group.add_argument(
        "--signal-min",
        type=float,
        help="Minimum value for pseudobulk track plot (default: data min)",
    )

    group.add_argument(
        "--signal-max",
        type=float,
        help="Maximum value for pseudobulk track plot (default: data max)",
    )

    group.add_argument(
        "--map-min",
        type=float,
        help="Minimum value for cell heatmap (default: data min)",
    )

    group.add_argument(
        "--map-max",
        type=float,
        help="Maximum value for cell heatmap (default: data max)",
    )

    group.add_argument(
        "--color",
        type=str,
        default="red",
        help="Color for top line plot (default: red)",
    )

    group.add_argument(
        "--colormap",
        type=str,
        default="Reds",
        help="Colormap for heatmap (default: Reds)",
    )

    group.add_argument(
        "--figsize",
        type=float,
        nargs=2,
        default=(14, 8),
        help="Figure size in inches (width height, default: 14 8)",
    )

    group.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="DPI for output PNG (default: 300)",
    )

    return parser


def parse_arguments(args=None):
    io_args = ParserCommon.inputOutputOptions(opts=["h5adfile", "outFile"], requiredOpts=["input", "outFile"])
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        prog="scPlotRegion",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[io_args, get_args(), other_args],
        description="""
``scPlotRegion`` plots pseudobulk chromatin signal across a genomic region with individual cell
profiles below.
""",
        usage="""
scPlotRegion -i INPUT_binned.h5ad --region CHROM:START-END -o OUTPUT_PREFIX
""",
        add_help=False,
    )

    # If no arguments are provided, show help and exit
    if args is None and len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args(args)

    # Validate mode is specified
    if args.mode is None:
        parser.error("You must specify --mode. Choose either 'activities' or 'aggregate'.")

    return args


def parse_region(region_str):
    """Parse region in either CHROM or CHROM:START-END format.

    Returns a tuple (chrom, start, end). For chromosome-only input `start` and
    `end` will be returned as 0 and np.inf so every feature in the chromosome is included.
    """
    s = region_str.strip()

    # CHROM only
    if re.fullmatch(r"[^:]+", s):
        return s, None, None

    # CHROM:START-END
    match = re.fullmatch(r"([^:]+):(\d+)-(\d+)", s)
    if not match:
        raise ValueError(f"Invalid region '{region_str}'. Expected 'CHROM' or 'CHROM:START-END'")

    chrom, start_str, end_str = match.groups()
    start = int(start_str)
    end = int(end_str)
    if start > end:
        raise ValueError(f"Invalid region '{region_str}'. START must be <= END.")

    return chrom, start, end


def get_region_mask(var_df, chrom: str, start: int, end: int):
    required_cols = {"chrom", "start", "end"}
    missing = required_cols.difference(var_df.columns)
    if missing:
        raise ValueError(f"anndata.var is missing required columns: {sorted(missing)}")

    # Keep features overlapping the query interval.
    return (
        (var_df["chrom"].astype(str) == str(chrom))
        & (var_df["end"].astype(int) >= start)
        & (var_df["start"].astype(int) <= end)
    )


# Helper function to reduce matrix width for heatmap plotting by summing contiguous columns.
def summarize_matrix_for_heatmap(X, dpi, figsize, heatmap_fraction=0.8):
    """
    Reduce matrix width for heatmap plotting by summing contiguous column bins.

    Parameters
    ----------
    X : array-like, shape (n_rows, n_cols)
        Input matrix.
    dpi : int or float
        Figure DPI used for plotting.
    figsize : tuple(float, float)
        Matplotlib figure size in inches.
    heatmap_fraction : float, default=0.8
        Fraction of horizontal figure width occupied by the heatmap.

    Returns
    -------
    X_binned : np.ndarray
        Matrix with same number of rows and reduced number of columns.
    bin_size : int
        Number of original columns summed per output column.
    max_cols : int
        Target maximum number of columns based on pixel budget.
    """
    X_arr = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    if X_arr.ndim != 2:
        raise ValueError("X must be a 2D matrix.")

    n_rows, n_cols = X_arr.shape
    max_cols = max(1, int(dpi * figsize[0] * heatmap_fraction))

    # No reduction needed
    if n_cols <= max_cols:
        return X_arr, 1, max_cols

    # Number of original bins merged into one plotted bin
    bin_size = int(np.ceil(n_cols / max_cols))
    n_new_cols = int(np.ceil(n_cols / bin_size))

    # Pad to a multiple of bin_size, then sum within each group
    pad_cols = n_new_cols * bin_size - n_cols
    if pad_cols > 0:
        X_arr = np.pad(X_arr, ((0, 0), (0, pad_cols)), mode="constant")

    X_binned = X_arr.reshape(n_rows, n_new_cols, bin_size).sum(axis=2)
    return X_binned, bin_size, max_cols


# Helper function to format the unit of coordinate tick labels.
def format_genomic_tick_labels(min_val: int, max_val: int, n_ticks: int = 9):
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
        tick_labels = [f"{v/scale:.{decimals}f} {unit}" for v in tick_vals]

    return tick_vals, tick_labels, unit


def main():
    args = get_args()

    chrom, region_start, region_end = parse_region(args.region)
    adata = ad.read_h5ad(args.input)

    # If region is CHROM-only (start/end are None) select all features on that
    # chromosome and infer bounds from the selected features. Otherwise use the
    # interval mask.
    if region_start is None or region_end is None:
        region_mask = adata.var["chrom"].astype(str) == str(chrom)
    else:
        region_mask = get_region_mask(adata.var, chrom, region_start, region_end)

    if region_mask.sum() == 0:
        raise ValueError(f"No features found in region {args.region}")

    var_region = adata.var.loc[region_mask].sort_values(["start", "end"])

    # If CHROM-only, infer start/end from features on that chromosome.
    if region_start is None or region_end is None:
        region_start = int(var_region["start"].min())
        region_end = int(var_region["end"].max())
    feature_idx = var_region.index
    adata_region = adata[:, feature_idx]

    X = adata_region.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    else:
        X = np.asarray(X)

    if X.ndim != 2:
        raise ValueError("Expected adata.X to be a 2D matrix.")

    if args.mode == "sum":
        agg = X.sum(axis=0)
    else:
        agg = X.mean(axis=0)

    x_positions = np.arange(X.shape[1])

    if args.summarize:
        X = summarize_matrix_for_heatmap(X, dpi=args.dpi, figsize=args.figsize)[0]

    # Reorder rows so higher-sum rows are on top
    row_order = np.argsort(X.sum(axis=1))[::-1]
    X_plot = X[row_order, :]

    fig, axes = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=args.figsize,
        sharex=True,
        gridspec_kw={"height_ratios": [1, 3]},
    )

    # Top line plot + fill
    axes[0].plot(x_positions, agg, color=args.color, linewidth=1.5)
    axes[0].fill_between(x_positions, agg, 0, color=args.color, alpha=0.9, linewidth=0)
    axes[0].set_ylabel(None)
    axes[0].tick_params(axis="x", which="both", bottom=False, labelbottom=False)

    if agg.min() == 0 and agg.max() == 0:
        axes[0].set_ylim(0, 1)

    # Only min/max ticks on top Y axis
    axes[0].set_ylim(agg.min(), agg.max())
    ymin, ymax = float(np.min(agg)), float(np.max(agg))
    axes[0].set_yticks([ymin, ymax])

    # Set Y limits based on user input or data min/max or arguments.
    if args.signal_min:
        signal_min = args.signal_min
    else:
        signal_min = float(agg.min())

    if args.signal_max:
        signal_max = args.signal_max
    else:
        signal_max = float(agg.max())

    # If all values are zero, set color limits to [0, 1].
    if X_plot.min == 0 and X_plot.max == 0:
        sys.output.write("All values in track plot are zero. Setting color limits to [0, 1].\n")
        map_min, map_max = 0, 1

    axes[0].set_ylim(signal_min, signal_max)

    # Set heatmap color limits based on user input or data min/max or arguments.
    if args.map_min:
        map_min = args.map_min
    else:
        map_min = float(X_plot.min())

    if args.map_max:
        map_max = args.map_max
    else:
        map_max = float(X_plot.max())

    # If all values are zero, set color limits to [0, 1].
    if X_plot.min == 0 and X_plot.max == 0:
        sys.output.write("All values in heatmap are zero. Setting color limits to [0, 1].\n")
        map_min, map_max = 0, 1

    # Heatmap
    im = axes[1].imshow(
        X_plot,
        aspect="auto",
        cmap=args.colormap,
        interpolation="nearest",
        vmin=map_min,
        vmax=map_max,
    )

    # Keep x-axis aligned with top plot
    axes[0].set_xlim(0, X.shape[1] - 1)

    n_ticks = min(9, int(args.figsize[1] / 2))
    n_ticks = min(n_ticks, X.shape[1])  # Avoid more ticks than features
    tick_positions = np.linspace(0, X.shape[1] - 1, n_ticks)
    _, tick_labels, unit = format_genomic_tick_labels(region_start, region_end, n_ticks=n_ticks)

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
    fig.tight_layout(rect=[0.08, 0, 1, 1])

    # Add colorbar to the LEFT of the heatmap
    hm_pos = axes[1].get_position()
    cax = fig.add_axes([hm_pos.x0 - 0.03, hm_pos.y0, 0.015, hm_pos.height])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(None)
    cbar.ax.yaxis.set_ticks_position("left")
    cbar.ax.yaxis.set_label_position("left")

    plt.margins(x=0.1, y=0.1)
    plt.show()
    fig.savefig(args.outputFilePrefix + ".png", dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
