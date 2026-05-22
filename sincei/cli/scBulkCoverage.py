#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import sys
import logging

from . import ParserCommon

# from sincei import _sincei

DESCRIPTION = """
``scBulkCoverage`` takes alignments of reads or fragments as input (BAM files), along with cell
grouping information, such as barcode -> batch, or barcode -> cluster, as tsv file, and generates a  coverage
track (bigWig or bedGraph) per group as output. The coverage is calculated as the number of reads per bin,
where bins are short consecutive counting windows of a defined  size. It is possible to extended/change the
length of the reads to better reflect the actual fragment length. ``scBulkCoverage`` offers normalization per
cluster using different methods.
"""

USAGE = "scBulkCoverage -b file1.bam file2.bam --labels file1 file2 -g group_info.tsv -o coverage"


def parse_arguments(args: list[str] | None = None) -> argparse.ArgumentParser:
    io_args = ParserCommon.inputOutputOptions(opts=["bamfiles", "groupInfo", "outFilePrefix"])
    bam_args = ParserCommon.bamOptions(default_opts={"binSize": 100})
    read_args = ParserCommon.readOptions()
    filter_args = ParserCommon.filterOptions()
    other_args = ParserCommon.otherOptions()
    parser = argparse.ArgumentParser(
        parents=[io_args, get_args(), bam_args, filter_args, read_args, other_args],
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


def get_args() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)

    optional = parser.add_argument_group("Coverage-related Options")
    optional.add_argument(
        "-of",
        "--outFileFormat",
        help='Output file type. Either "bigwig" or "bedgraph".',
        choices=["bigwig", "bedgraph"],
        default="bigwig",
    )

    optional.add_argument(
        "-n",
        "--normalizeUsing",
        help="How to normalize the pseudo-bulk counts. Options are "
        '"CPM": normalized each bin to the counts per million mapped reads in that group.\n'
        '"Frequency": binarize the coverage per bin and normalize to the total no. of cells per group. \n'
        '"Mean": get mean signal per bin across cells in each group.\n'
        '"None": simply return the sum of coverage per group.',
        choices=["CPM", "Frequency", "Mean", "None"],
        default="CPM",
    )

    optional.add_argument(
        "-ig",
        "--ignoreForNormalization",
        metavar="CHR",
        help="Chromosomes to skip while calculating normalization factors",
        nargs="+",
        default=None,
    )

    optional.add_argument(
        "-nr",
        "--normalizeByReference",
        help="NOT IMPLEMENTED: Normalize each group of cells by a reference group (which must be present in the --groupinfo file)"
        "Note that the --normalizeUsing method is applied beforehand.",
        choices=["ratio", "log2_ratio", "difference", "None"],
        default=None,
    )

    optional.add_argument(
        "--scaleFactor",
        metavar="FLOAT",
        help="The computed scaling factor (or 1, if not applicable) will "
        "be multiplied by this. (Default: %(default)s)",
        default=1.0,
        type=float,
        required=False,
    )

    optional.add_argument(
        "--MNase",
        help="Determine nucleosome positions from MNase-seq/CUTnRUN data. "
        "Only 3 nucleotides at the center of each fragment are counted. "
        "The fragment ends are defined by the two mate reads. Only fragment lengths"
        "between 130 - 200 bp are considered to avoid dinucleosomes or other artifacts. "
        "By default, any fragments smaller or larger than this are ignored. To "
        "over-ride this, use the --minFragmentLength and --maxFragmentLength options, "
        "which will default to 130 and 200 if not otherwise specified in the presence "
        "of --MNase. *NOTE*: Requires paired-end data. A bin size of 1 is recommended.",
        action="store_true",
    )

    optional.add_argument(
        "--Offset",
        metavar="INT",
        help="Uses this offset inside of each read as the signal. This is useful in "
        "cases like RiboSeq or GROseq, where the signal is 12, 15 or 0 bases past the "
        "start of the read. This can be paired with the --filterRNAstrand option. "
        "Note that negative values indicate offsets from the end of each read. A value "
        "of 1 indicates the first base of the alignment (taking alignment orientation "
        "into account). Likewise, a value of -1 is the last base of the alignment. An "
        "offset of 0 is not permitted. If two values are specified, then they will be "
        "used to specify a range of positions. Note that specifying something like "
        "--Offset 5 -1 will result in the 5th through last position being used, which "
        "is equivalent to trimming 4 bases from the 5-prime end of alignments. Note "
        "that if you specify --centerReads, the centering will be performed before the "
        "offset.",
        type=int,
        nargs="+",
        required=False,
    )

    return parser


def main(args: list[str] | None = None) -> int:
    args = parse_arguments().parse_args(args)

    for arg in args.__dict__:
        logging.debug(f"{arg}: {args.__dict__[arg]}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
