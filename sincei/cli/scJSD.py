#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import sys
import logging

from . import ParserCommon

# from sincei import _sincei

DESCRIPTION = """
``scJSD`` samples regions in the genome from BAM files and compares the cumulative read coverages
for each cell on those regions to a synthetic cell with poisson distributed reads using the Jensen-Shannon
distance. Cells with high enrichment of signals show a higher JSD compared to cells whose signal is
homogeneously distributed.
"""

USAGE = "scJSD -b treatment.bam control.bam -bc barcodes.txt -o scJSD.tsv"


def parse_arguments(args: list[str] | None = None) -> argparse.ArgumentParser:
    io_args = ParserCommon.inputOutputOptions(
        opts=["bamfiles", "barcodes", "outFile"],
        requiredOpts=["bamfiles", "barcodes", "outFile"],
    )
    bam_args = ParserCommon.bamOptions(suppress_args=["region", "distanceBetweenBins"], default_opts={"binSize": 10000})
    read_args = ParserCommon.readOptions(suppress_args=["filterRNAstrand", "extendReads", "centerReads"])
    filter_args = ParserCommon.filterOptions(
        suppress_args=[
            "motifFilter",
            "genome2bit",
            "GCcontentFilter",
            "minAlignedFraction",
        ]
    )
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        parents=[io_args, bam_args, read_args, filter_args, get_args(), other_args],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=DESCRIPTION,
        conflict_handler="resolve",
        usage=USAGE,
        add_help=False,
    )

    # If no arguments are provided, show help and exit
    if args is None and len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    return parser


def get_args() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False, conflict_handler="resolve")
    optional = parser.add_argument_group("Optional arguments")

    optional.add_argument(
        "-n",
        "--numberOfSamples",
        metavar="INT",
        help="The number of bins that are sampled from the genome, "
        "for which the overlapping number of reads is computed. (Default: %(default)s)",
        default=1e5,
        type=int,
    )

    optional.add_argument(
        "--skipZeros",
        help="If set, then regions with zero overlapping reads"
        "for *all* given BAM files are ignored. This "
        "will result in a reduced number of read "
        "counts than that specified in --numberOfSamples",
        action="store_true",
    )

    return parser


def main(args: list[str] | None = None) -> int:
    args = parse_arguments().parse_args(args)

    for arg in args.__dict__:
        logging.debug(f"{arg}: {args.__dict__[arg]}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
