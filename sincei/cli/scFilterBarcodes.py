#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import sys
import logging

from . import ParserCommon

# from sincei import _sincei

DESCRIPTION = """
``scFilterBarcodes`` identifies barcodes present in a BAM file and produces a list. You can
optionally filter these barcodes by matching them to a whitelist or based on total counts.
"""

USAGE = "scFilterBarcodes -b sample.bam -w whitelist.txt -o barcodes_detected.txt"


def parse_arguments(args: list[str] | None = None) -> argparse.ArgumentParser:
    io_args = ParserCommon.inputOutputOptions(opts=["bamfile", "whitelist", "outFile"], requiredOpts=["bamfile"])
    bam_args = ParserCommon.bamOptions(
        suppress_args=["labels", "smartLabels", "distanceBetweenBins", "region"],
        default_opts={"binSize": 100000},
    )
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        parents=[io_args, get_args(), bam_args, other_args],
        formatter_class=argparse.RawDescriptionHelpFormatter,
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
    general = parser.add_argument_group("Counting Options")

    general.add_argument(
        "-d",
        "--minHammingDist",
        help="Minimum hamming distance to match the barcode in whitelist. Note that increasing the "
        "hamming distance really slows down the barcode detection process.",
        metavar="INT",
        type=int,
        default=0,
        required=False,
    )

    general.add_argument(
        "-mc",
        "--minCount",
        help="Minimum no. of bins with non-zero counts, in order to report a barcode. Note that this number would range "
        "from 0 to genome size/binSize. ",
        metavar="INT",
        type=int,
        default=0,
        required=False,
    )

    general.add_argument(
        "-mq",
        "--minMappingQuality",
        metavar="INT",
        type=int,
        help="If set, only reads that have a mapping quality score of at least this are considered.",
    )

    general.add_argument(
        "-rp",
        "--rankPlot",
        metavar="STR",
        type=str,
        help='The output file name to plot the ranked counts per barcode (similar to the "knee plot",'
        "but counts are be the number of non-zero bins in this case).",
    )

    return parser


def main(args: list[str] | None = None) -> int:
    args = parse_arguments().parse_args(args)

    for arg in args.__dict__:
        logging.debug(f"{arg}: {args.__dict__[arg]}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
