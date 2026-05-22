#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import sys
import logging

from . import ParserCommon

DESCRIPTION = """
``scExportSignal`` export AnnData to other formats.
"""

USAGE = "scExportSignal -i INPUT.h5ad -f FORMAT -o OUTPUT"


def parse_arguments(args: list[str] | None = None) -> argparse.ArgumentParser:
    io_args = ParserCommon.inputOutputOptions(
        opts=["h5adfile", "outFile", "region"], requiredOpts=["input", "outFilePrefix", "region"]
    )
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        prog="scExportSignal",
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

    args = parser.parse_args(args)

    return args


def get_args() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)

    group = parser.add_argument_group("Export Options")

    group.add_argument(
        "--outFileFormat",
        type=str,
        default="h5ad",
        choices=["bm", "mtx", "loom"],
        help="Output file format for `scExportSignal`. "
        '"bm" refers to the BedGraphMatrix format, useful for single-cell data visualization '
        "with pyGenomeTracks. It stores data in dense format, so we recommend choosing a "
        "region to export instead of the whole dataset. "
        '"mtx" refers to the MatrixMarket sparse-matrix format. The output in this case would '
        "be <prefix>.counts.mtx, along with <prefix>.rownames.txt and <prefix>.colnames.txt"
        '"loom" refers to the loom file format, an hddf5-based legacy format for single-cell '
        "genomics data.",
        required=True,
    )

    return parser


def main(args: list[str] | None = None) -> int:
    args = parse_arguments().parse_args(args)

    for arg in args.__dict__:
        logging.debug(f"{arg}: {args.__dict__[arg]}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
