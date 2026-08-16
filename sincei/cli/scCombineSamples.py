#!/usr/bin/env python
from __future__ import annotations

import sys
from typing import TYPE_CHECKING

from . import ParserCommon
from . import _parsers as backend

if TYPE_CHECKING:
    import argparse

DESCRIPTION = (
    "Concatenate/merge AnnData files from different "
    "samples.\n\n``scCombineSamples`` combines multiple count matrices "
    "(output of scCountReads) into one. Only features present in all matrices "
    "will be kept. The result is a .h5ad (AnnData) file containing the "
    "combined count matrix.\n*NOTE*: it doesn't perform any 'batch effect "
    "correction' or 'integration' of data from different technologies."
)

USAGE = "scCombineSamples -i a.h5ad b.h5ad -o combined.h5ad"


def parse_arguments() -> argparse.ArgumentParser:
    return ParserCommon.build_parser(
        DESCRIPTION,
        USAGE,
        [
            ParserCommon.input_output_options(["h5ad_files", "out_file"]),
            ParserCommon.bam_options(
                ["labels", "smart_labels"], panel="Input / Output options"
            ),
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
        labels=args.labels,
        out_file=args.outFile,
        smart_labels=args.smartLabels,
        number_of_processors=args.numberOfProcessors,
        verbose=args.verbose,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
