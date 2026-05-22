#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import sys
import logging

from . import ParserCommon

DESCRIPTION = """
``scPlotRegion`` plots pseudobulk chromatin signal across a genomic region with individual cell
profiles below.
"""

USAGE = "scPlotRegion -i binned_data.h5ad --region CHROM:START-END -o output.png"


def parse_arguments(args: list[str] | None = None) -> argparse.ArgumentParser:
    io_args = ParserCommon.inputOutputOptions(
        opts=["h5adfile", "outFile", "region"], requiredOpts=["input", "outFile", "region"]
    )
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        prog="scPlotRegion",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[io_args, get_args(), other_args],
        description=DESCRIPTION,
        usage=USAGE,
        add_help=False,
    )

    # If no arguments are provided, show help and exit
    if args is None and len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    return parser


def get_args() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group("Plot Options")

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
        metavar="FLOAT",
        type=float,
        help="Minimum value for pseudobulk track plot (default: data min)",
    )

    group.add_argument(
        "--signal-max",
        metavar="FLOAT",
        type=float,
        help="Maximum value for pseudobulk track plot (default: data max)",
    )

    group.add_argument(
        "--map-min",
        metavar="FLOAT",
        type=float,
        help="Minimum value for cell heatmap (default: data min)",
    )

    group.add_argument(
        "--map-max",
        metavar="FLOAT",
        type=float,
        help="Maximum value for cell heatmap (default: data max)",
    )

    group.add_argument(
        "--color",
        metavar="STR",
        type=str,
        default="red",
        help="Color for top line plot (default: red)",
    )

    group.add_argument(
        "--colormap",
        metavar="STR",
        type=str,
        default="Reds",
        help="Colormap for heatmap (default: Redss)",
    )

    group.add_argument(
        "--figsize",
        metavar="INT",
        type=float,
        nargs=2,
        default=(14, 8),
        help="Figure size in inches (width height, default: 14 8)",
    )

    group.add_argument(
        "--dpi",
        metavar="INT",
        type=int,
        default=300,
        help="DPI for output PNG (default: 300)",
    )

    return parser


def main(args: list[str] | None = None) -> int:
    args = parse_arguments().parse_args(args)

    for arg in args.__dict__:
        logging.debug(f"{arg}: {args.__dict__[arg]}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
