#!/usr/bin/env python
from __future__ import annotations

import sys
from typing import TYPE_CHECKING

from . import ParserCommon
from . import _parsers as backend

if TYPE_CHECKING:
    import argparse

DESCRIPTION = (
    "Merge AnnData files from different modalities into a MuData "
    "object.\n\n``scCombineMods`` combines multiple count matrices (output of "
    "scCountReads) of different data modalities (e.g. gene expression, "
    "chromatin accessibility, histone modifications) into one. The result is "
    "a .h5mu (MuData) file containing each of the data modalities "
    "provided.\n*NOTE*: this does not perform any 'batch effect correction' "
    "or 'integration' of data from different technologies."
)

USAGE = "scCombineMods -i rna.h5ad atac.h5ad -l RNA ATAC -o combined.h5mu"


def parse_arguments() -> argparse.ArgumentParser:
    return ParserCommon.build_parser(
        DESCRIPTION,
        USAGE,
        [
            ParserCommon.input_output_options(["h5ad_files", "out_file"]),
            ParserCommon.bam_options(["labels"], defaults={"labels_required": True}),
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
        number_of_processors=args.numberOfProcessors,
        verbose=args.verbose,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
