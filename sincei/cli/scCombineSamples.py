#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import sys
import logging

from . import ParserCommon

DESCRIPTION = """
``scCombineSamples`` combines multiple count matrices (output of scCountReads) into one. Only
features present in all matrices will be kept. The result is a .h5ad (AnnData) file containing the
combined count matrix.
*NOTE*: it doesn't perform any 'batch effect correction' or 'integration' of data from different
technologies.
"""

USAGE = "scCombineSamples -i sample1.h5ad sample2.h5ad -o combined.h5ad"


def parse_arguments(args: list[str] | None = None) -> argparse.ArgumentParser:
    io_args = ParserCommon.inputOutputOptions(
        opts=["h5adfiles", "labels", "outFile"], requiredOpts=["h5adfiles", "labels", "outFile"]
    )
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        parents=[io_args, other_args],
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


def main(args: list[str] | None = None) -> int:
    args = parse_arguments().parse_args(args)

    for arg in args.__dict__:
        logging.debug(f"{arg}: {args.__dict__[arg]}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
