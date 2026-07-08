#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re
import sys

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np

from sincei import ParserCommon

### ------ Functions ------


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


### ------ Args ------


def parseArguments(args=None):
    io_args = ParserCommon.inputOutputOptions(
        opts=["h5adfile", "outFilePrefix"], requiredOpts=["input", "outFilePrefix"]
    )
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        prog="scPlotRegion",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[io_args, get_args(), other_args],
        description="""
``scPlotRegion`` plots signal from individual cells in a given genomic region as a heatmap (in the style of a trackplot),
togather with the summary profiles (pseudobulk signal) on top.
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

    return parser


def get_args():
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group("Plot Options")

    group.add_argument(
        "--region",
        "-r",
        required=True,
        help="Genomic region in CHROM:START-END format",
    )

    group.add_argument(
        "--mode",
        "-m",
        required=False,
        choices=["sum", "mean"],
        default="sum",
        help="How to aggregate signal for plotting the summary profiles (pseudobulk) on top. (default: %(default)s)",
    )

    group.add_argument(
        "--signalMin",
        type=float,
        help="Minimum value for pseudobulk track plot (default: minimum signal in the region)",
    )

    group.add_argument(
        "--signalMax",
        type=float,
        help="Maximum value for pseudobulk summary profile (default: maximum signal in the region)",
    )

    group.add_argument(
        "--hmapMin",
        type=float,
        help="Minimum value for single-cell heatmap signal (default: minimum signal in the region)",
    )

    group.add_argument(
        "--hmapMax",
        type=float,
        help="Maximum value for single-cell heatmap signal (default: maximum signal in the region)",
    )

    group.add_argument(
        "--summaryPlotColor",
        type=str,
        default="red",
        help="Color for the pseudobulk summary profile plot (default: %(default)s)",
    )

    group.add_argument(
        "--colorMap",
        type=str,
        default="Reds",
        help="Colormap for heatmap. Must be a valid `matplotlib colormap` entry (default: %(default)s)",
    )

    group.add_argument(
        "--figsize",
        type=float,
        nargs=2,
        default=(14, 8),
        help="Figure size in inches (<width>, <height>; (default: %(default)s))",
    )

    group.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="DPI for output PNG (default: %(default)s)",
    )

    group.add_argument(
        "--outFileFormat",
        help="image format type. If given, this option "
        "overrides the image format inferred from the suffix of the outFilePrefix",
        default=None,
        choices=["png", "pdf", "svg", "eps", "plotly"],
    )

    return parser


def main(args=None):
    args = parseArguments().parse_args(args)

    chrom, region_start, region_end = ParserCommon.parse_region(args.region)
    adata = ParserCommon.validateAnndata(ad.read_h5ad(args.input), args.input)

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
    axes[0].plot(x_positions, agg, color=args.summaryPlotColor, linewidth=1.5)
    axes[0].fill_between(x_positions, agg, 0, color=args.summaryPlotColor, alpha=0.9, linewidth=0)
    axes[0].set_ylabel(None)
    axes[0].tick_params(axis="x", which="both", bottom=False, labelbottom=False)

    ymin, ymax = float(np.min(agg)), float(np.max(agg))
    if ymin == 0.0 and ymax == 0.0:
        ymax = 1.0

    # Set Y limits based on user input or data min/max.
    if args.signalMin is not None:
        ymin = args.signalMin
    if args.signalMax is not None:
        ymax = args.signalMax

    axes[0].set_ylim(ymin, ymax)
    # Only min/max ticks on top Y axis
    axes[0].set_yticks([ymin, ymax])
    # Set heatmap color limits based on user input or data min/max.
    if args.hmapMin is not None:
        map_min = args.hmapMin
    else:
        map_min = float(X_plot.min())

    if args.hmapMax is not None:
        map_max = args.hmapMax
    else:
        map_max = float(X_plot.max())

    # If all values are zero, set color limits to [0, 1] so imshow has a valid range.
    if map_min == 0.0 and map_max == 0.0:
        sys.stderr.write("All values in heatmap are zero. Setting color limits to [0, 1].\n")
        map_max = 1.0

    # Heatmap
    im = axes[1].imshow(
        X_plot,
        aspect="auto",
        cmap=args.colorMap,
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

    if not args.outFileFormat:
        args.outFileFormat = None
    plt.margins(x=0.1, y=0.1)
    fig.savefig(args.outFilePrefix, dpi=args.dpi, bbox_inches="tight", format=args.outFileFormat)
    plt.close(fig)


if __name__ == "__main__":
    main()
