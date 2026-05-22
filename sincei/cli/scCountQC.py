#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import sys
import logging

from . import ParserCommon

# from sincei import _sincei

DESCRIPTTION = """
``scCountQC`` calculates multiple quality controls metrics on the input .h5ad file (output of scCountReads) and
(optionally) filters the input file based on filterCellArgs/filterRegionArgs. The output is either an updated .h5ad
object (if filtering is requested) or the filtering metrics (--outMetrics).
"""

USAGE = "scCountQC -i cellCounts.h5ad -o cellCounts.filtered.h5ad -om qc_metrics.tsv"


def parse_arguments(args: list[str] | None = None) -> argparse.ArgumentParser:
    io_args = ParserCommon.inputOutputOptions(opts=["h5adfile", "outFile"])
    other_args = ParserCommon.otherOptions()
    parser = argparse.ArgumentParser(
        parents=[io_args, get_args(), other_args],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=DESCRIPTTION,
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

    general = parser.add_argument_group("Filtering arguments")
    general.add_argument(
        "-d",
        "--describe",
        action="store_true",
        help="Print a list of cell and region metrics available for QC/filtering.",
    )

    general.add_argument(
        "-om",
        "--outMetrics",
        type=str,
        help="Prefix of the output file with calculated QC metrics. If given, the cell metrics are printed in "
        "<prefix>.cell.tsv and region metrics as <prefix>.region.tsv",
    )

    general.add_argument(
        "-fc",
        "--filterCellArgs",
        type=str,
        help='List of arguments to filter cells. The format is "arg_name: minvalue, maxvalue; arg_name: minvalue, maxvalue; ...." '
        "where arg_name is the QC metric for cells present in the input h5ad file. In order to view all available "
        'cell filtering metrics, run scCountFilter with "--describe". The two arguments are supplied (minvalue, maxvalue) '
        "they are used as lower and upper bounds to filter cells. Make sure they are float/integer numbers.",
    )

    general.add_argument(
        "-fr",
        "--filterRegionArgs",
        type=str,
        help='List of arguments to filter regions. The format is "arg_name: minvalue, maxvalue; arg_name: minvalue; ...." '
        "where arg_name is the QC metric for regions present in the input h5ad file. In order to view all available "
        'cell filtering metrics, run scCountFilter with "--describe". The two arguments are supplied (minvalue, maxvalue) '
        "they are used as lower and upper bounds to filter cells. Make sure they are float/integer numbers.",
    )

    general.add_argument(
        "-rb",
        "--region_blacklist",
        help="A BED or GTF file containing regions that should be excluded from all analyses. "
        "Regions in the anndata object that overlap with blacklisted regions will be removed.",
        metavar="BED",
        nargs="+",
    )

    general.add_argument(
        "-cb",
        "--cell_blacklist",
        default=None,
        type=argparse.FileType("r"),
        help="A list of barcodes to be excluded for the clustering. The barcodes "
        "(along with sample labels) must be present in the input object.",
    )

    general.add_argument(
        "-chb",
        "--chrom_blacklist",
        default=None,
        help="A space separated list of chromosomes to exclude. eg. chrM chrUn",
        metavar="CHR1",
        nargs="+",
    )

    return parser


def main(args: list[str] | None = None) -> int:
    args = parse_arguments().parse_args(args)

    for arg in args.__dict__:
        logging.debug(f"{arg}: {args.__dict__[arg]}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
