#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import sys
import logging

from . import ParserCommon

# from sincei import _sincei


DESCRIPTION = """
``scFilterStats`` estimates the number of reads that would be filtered given a set of criteria
and prints it to the terminal. Furthermore, it tracks the number of singleton reads.
The following metrics will always be tracked regardless of what you specify (the order output also matches this):

* Total reads (including unmapped)
* Mapped reads
* Reads in blacklisted regions (--blackListFileName)

The following metrics are estimated according to the --binSize and --distanceBetweenBins parameters:

* Estimated mapped reads filtered (the total number of mapped reads filtered for any reason)
* Alignments with a below threshold MAPQ (--minMappingQuality)
* Alignments with at least one missing flag (--samFlagInclude)
* Alignments with undesirable flags (--samFlagExclude)
* Duplicates determined by sincei (--duplicateFilter)
* Duplicates marked externally (e.g., by picard)
* Singletons (paired-end reads with only one mate aligning)
* Wrong strand (due to --filterRNAstrand)

The sum of these may be more than the total number of reads. Note that alignments are sampled from
bins of size --binSize spaced --distanceBetweenBins apart.
"""

USAGE = "scFilterStats -b sample1.bam sample2.bam -bc barcodes.txt -bl blacklist.bed -o stats.tsv"


def parse_arguments(args: list[str] | None = None) -> argparse.ArgumentParser:
    io_args = ParserCommon.inputOutputOptions(opts=["bamfiles", "barcodes", "outFile"], requiredOpts=["barcodes"])
    bam_args = ParserCommon.bamOptions(
        suppress_args=["region", "groupTag"],
        default_opts={"binSize": 100000, "distanceBetweenBins": 1000000},
    )
    filter_args = ParserCommon.filterOptions()
    read_args = ParserCommon.readOptions(
        suppress_args=[
            "minFragmentLength",
            "maxFragmentLength",
            "extendReads",
            "centerReads",
        ]
    )
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        parents=[io_args, bam_args, filter_args, read_args, other_args],
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
