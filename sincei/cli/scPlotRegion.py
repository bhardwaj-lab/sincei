#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys

from . import ParserCommon
from . import _parsers as backend

DESCRIPTION = "Plot pseudo-bulk and per cell coverage for a genomic region."

USAGE = "scPlotRegion -i counts.h5ad -r chr1:1000:2000 -m mean -o plot.png"


_DISPLAY = "Display options"
_COLOR = "Color / Scale options"
_FIGURE = "Figure options"


def display_options() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    display = parser.add_argument_group(_DISPLAY)
    display.add_argument(
        "-m",
        "--mode",
        required=True,
        help="Aggregation mode for the top subplot.",
    )
    color = parser.add_argument_group(_COLOR)
    color.add_argument(
        "--signalMin",
        dest="signalMin",
        type=float,
        default=None,
        help="Minimum value for pseudobulk track plot (default: data min).",
    )
    color.add_argument(
        "--signalMax",
        dest="signalMax",
        type=float,
        default=None,
        help="Maximum value for pseudobulk track plot (default: data max).",
    )
    color.add_argument(
        "--mapMin",
        dest="mapMin",
        type=float,
        default=None,
        help="Minimum value for cell heatmap (default: data min).",
    )
    color.add_argument(
        "--mapMax",
        dest="mapMax",
        type=float,
        default=None,
        help="Maximum value for cell heatmap (default: data max).",
    )
    color.add_argument(
        "--color",
        default="red",
        help="Color for the top line plot.",
    )
    color.add_argument(
        "--colormap",
        default="Reds",
        help="Colormap for the heatmap.",
    )
    figure = parser.add_argument_group(_FIGURE)
    figure.add_argument(
        "--figWidth",
        dest="figWidth",
        type=float,
        default=14.0,
        help="Figure width in inches.",
    )
    figure.add_argument(
        "--figHeight",
        dest="figHeight",
        type=float,
        default=8.0,
        help="Figure height in inches.",
    )
    figure.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="DPI for the output PNG.",
    )
    return parser


def parse_arguments() -> argparse.ArgumentParser:
    return ParserCommon.build_parser(
        DESCRIPTION,
        USAGE,
        [
            ParserCommon.input_output_options(
                ["h5ad_file", "out_file", "region"],
                defaults={"region_required": True},
            ),
            display_options(),
            ParserCommon.other_options(),
        ],
    )


def main(argv: list[str] | None = None) -> int:
    parser = parse_arguments()
    if argv is None and len(sys.argv) == 1:
        parser.print_help()
        return 0
    args = parser.parse_args(argv)

    backend.log_parameters(
        input=args.input,
        out_file=args.outFile,
        region=args.region,
        mode=args.mode,
        signal_min=args.signalMin,
        signal_max=args.signalMax,
        map_min=args.mapMin,
        map_max=args.mapMax,
        color=args.color,
        colormap=args.colormap,
        fig_width=args.figWidth,
        fig_height=args.figHeight,
        dpi=args.dpi,
        number_of_processors=args.numberOfProcessors,
        verbose=args.verbose,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
