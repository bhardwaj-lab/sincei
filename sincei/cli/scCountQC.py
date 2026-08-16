#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys

from . import ParserCommon
from . import _parsers as backend

DESCRIPTION = (
    "Perform quality control and filter cells and regions from a "
    "cell-by-feature matrix.\n\n``scCountQC`` calculates multiple quality "
    "controls metrics on the input .h5ad file (output of scCountReads) and "
    "(optionally) filters the input file based on "
    "filterCellArgs/filterRegionArgs. The output is either an updated .h5ad "
    "object (if filtering is requested) or the filtering metrics "
    "(--outMetrics)."
)

USAGE = "scCountQC -i counts.h5ad -o filtered.h5ad"


_QC = "QC options"


def qc_options() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(_QC)
    group.add_argument(
        "-d",
        "--describe",
        action="store_true",
        help="Print a list of cell and region metrics available for QC/filtering.",
    )
    group.add_argument(
        "-om",
        "--outMetrics",
        dest="outMetrics",
        default=None,
        help="Prefix of the output file with calculated QC metrics. If given, the "
        "cell metrics are printed in <prefix>.cell.tsv and region metrics as "
        "<prefix>.region.tsv.",
    )
    group.add_argument(
        "-fc",
        "--filterCellArgs",
        dest="filterCellArgs",
        default=None,
        help='List of arguments to filter cells. The format is "arg_name: minvalue, '
        'maxvalue; arg_name: minvalue, maxvalue; ...." where arg_name is a cell QC '
        'metric present in the input h5ad file. Run with "--describe" to view all '
        "available metrics. The two values are used as lower and upper bounds to "
        "filter cells.",
    )
    group.add_argument(
        "-fr",
        "--filterRegionArgs",
        dest="filterRegionArgs",
        default=None,
        help='List of arguments to filter regions. The format is "arg_name: minvalue, '
        'maxvalue; arg_name: minvalue; ...." where arg_name is a region QC metric '
        'present in the input h5ad file. Run with "--describe" to view all available '
        "metrics. The two values are used as lower and upper bounds to filter regions.",
    )
    group.add_argument(
        "-rb",
        "--regionBlacklist",
        dest="regionBlacklist",
        metavar="BED",
        action="extend",
        nargs="+",
        default=None,
        help="A BED or GTF file containing regions that should be excluded from all "
        "analyses. Regions in the anndata object that overlap with blacklisted regions "
        "will be removed.",
    )
    group.add_argument(
        "-cb",
        "--cellBlacklist",
        dest="cellBlacklist",
        default=None,
        help="A list of barcodes to be excluded from the clustering. The barcodes "
        "(along with sample labels) must be present in the input object.",
    )
    group.add_argument(
        "-chb",
        "--chromBlacklist",
        dest="chromBlacklist",
        metavar="CHR",
        action="extend",
        nargs="+",
        default=None,
        help="A space separated list of chromosomes to exclude, e.g. chrM chrUn.",
    )
    return parser


def parse_arguments() -> argparse.ArgumentParser:
    return ParserCommon.build_parser(
        DESCRIPTION,
        USAGE,
        [
            ParserCommon.input_output_options(["h5ad_file", "out_file"]),
            qc_options(),
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
        describe=args.describe,
        out_metrics=args.outMetrics,
        filter_cell_args=args.filterCellArgs,
        filter_region_args=args.filterRegionArgs,
        region_blacklist=args.regionBlacklist,
        cell_blacklist=args.cellBlacklist,
        chrom_blacklist=args.chromBlacklist,
        number_of_processors=args.numberOfProcessors,
        verbose=args.verbose,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
