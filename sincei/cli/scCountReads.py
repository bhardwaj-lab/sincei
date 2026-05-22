#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import sys
import logging

from . import ParserCommon

# from sincei import _sincei

DESCRIPTION = """
``scCountReads`` computes the read coverages per cell barcode for genomic regions in the provided
BAM file(s). The analysis can be performed for the entire genome by running the program in ``bins`` mode. If
you want to count the read coverage for specific regions only, use the ``features`` mode instead. The
standard output of ``scCountReads`` is a ".h5ad" file with counts, along with rowName (features) and colNames
(cell barcodes).

Detailed help for each sub-command is available by typing:

    scCountReads bins -h
    scCountReads features -h
"""

USAGE = "scCountReads <method> --bamfiles file1.bam file2.bam --barcodes whitelist.txt -o counts.h5ad [opts]"

METHODS = ["bins", "features"]


def parse_arguments(args: list[str] | None = None) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=DESCRIPTION,
        usage=USAGE,
        conflict_handler="resolve",
    )

    subparsers = parser.add_subparsers(title="Counting mode")

    read_args = ParserCommon.readOptions(suppress_args=["filterRNAstrand"])
    filter_args = ParserCommon.filterOptions()
    other_args = ParserCommon.otherOptions()

    # bins mode options
    subparsers.add_parser(
        "bins",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[
            ParserCommon.inputOutputOptions(
                opts=["bamfiles", "barcodes", "outFile", "BED"],
                requiredOpts=["barcodes", "outFile"],
                suppress_args=["BED"],
            ),
            ParserCommon.gtf_options(),
            ParserCommon.bamOptions(default_opts={"binSize": 10000, "distanceBetweenBins": 0}),
            read_args,
            filter_args,
            get_args(),
            other_args,
        ],
        help="The reads are counted in bins of equal size. The bin size and distance between bins can be adjusted.",
        add_help=False,
        usage="scCountReads bins -bs 10000 --bamfiles file1.bam file2.bam --barcodes whitelist.txt -o counts.h5ad",
    )

    # BED file arguments
    subparsers.add_parser(
        "features",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[
            ParserCommon.inputOutputOptions(
                opts=["bamfiles", "barcodes", "outFile", "BED"],
                requiredOpts=["barcodes", "outFile", "BED"],
            ),
            ParserCommon.gtf_options(),
            ParserCommon.bamOptions(
                suppress_args=["binSize", "distanceBetweenBins"],
                default_opts={"binSize": 10000, "distanceBetweenBins": 0},
            ),
            read_args,
            filter_args,
            get_args(),
            other_args,
        ],
        help="The user provides a BED/GTF file containing all regions that "
        "should be counted. A common use would be to count scRNA-seq reads on Genes.",
        usage="scCountReads features --BED selection.bed [genes.gtf] --bamfiles file1.bam file2.bam --barcodes whitelist.txt -o counts.h5ad",
        add_help=False,
    )

    # If no arguments are provided, show help and exit
    if args is None and len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    # If only subcommand is provided, show subcommand help
    if sys.argv[1] in METHODS and len(sys.argv) == 2:
        parser.parse_args([sys.argv[1], "-h"])

    return parser


def get_args() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    optional = parser.add_argument_group("Misc arguments")

    optional.add_argument(
        "--valueTag",
        metavar="STR",
        type=str,
        default=None,
        help='Instead of counting each read/fragment as "1", add the values '
        "from a given BAM tag to the count matrix. For example, this can be "
        "used to count the number of methylated CpG per fragment.",
        nargs="?",
    )

    optional.add_argument(
        "--genomeChunkSize",
        metavar="INT",
        type=int,
        default=None,
        help="Manually specify the size of the genome provided to each processor. "
        "The default value of None specifies that this is determined by read "
        "density of the BAM file.",
    )

    return parser


def main(args: list[str] | None = None) -> int:
    args = parse_arguments().parse_args(args)

    for arg in args.__dict__:
        logging.debug(f"{arg}: {args.__dict__[arg]}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
