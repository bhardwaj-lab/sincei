#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import sys
import logging

from . import ParserCommon

DESCRIPTION = """
``scCombineMods`` combines multiple count matrices (output of scCountReads) of different data
modalities (e.g. gene expression, chromatin accessibility, histone modifications) into one.
The result is a .h5mu (MuData) file containing each of the data modalities provided.
*NOTE*: it doesn't perform any 'batch effect correction' or 'integration' of data from different
technologies, which requires more sophisticated methods.
"""

USAGE = "scCombineMods -i modality1.h5ad modality2.h5ad -o combined.h5mu"


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
