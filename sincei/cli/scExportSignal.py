#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys

from . import ParserCommon
from . import _parsers as backend

DESCRIPTION = "Export .h5ad objects to other formats."

USAGE = "scExportSignal -i counts.h5ad --outFileFormat mtx -o export"


_EXPORT = "Export options"


def export_options() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(_EXPORT)
    group.add_argument(
        "--outFileFormat",
        dest="outFileFormat",
        metavar="FORMAT",
        choices=ParserCommon.EXPORT_FORMATS,
        required=True,
        help="Output file format.\n\n"
        "bm: BedGraphMatrix format, useful for single-cell visualization with "
        "pyGenomeTracks; stores data densely, so prefer exporting a region rather than "
        "the whole dataset.\n\n"
        "mtx: MatrixMarket sparse format (<prefix>.counts.mtx plus "
        "<prefix>.rownames.txt and <prefix>.colnames.txt).\n\n"
        "loom: loom hdf5-based legacy format for single-cell genomics data.",
    )
    return parser


def parse_arguments() -> argparse.ArgumentParser:
    return ParserCommon.build_parser(
        DESCRIPTION,
        USAGE,
        [
            ParserCommon.input_output_options(["h5ad_file", "out_prefix", "region"]),
            export_options(),
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
        out_file=args.outFilePrefix,
        out_file_format=args.outFileFormat,
        region=args.region,
        number_of_processors=args.numberOfProcessors,
        verbose=args.verbose,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
