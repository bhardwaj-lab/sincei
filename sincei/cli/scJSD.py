#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys

from . import ParserCommon
from . import _parsers as backend

DESCRIPTION = (
    "Compare read coverages on sampled regions using the Jensen-Shannon "
    "distance.\n\n``scJSD`` samples regions in the genome from BAM files and "
    "compares the cumulative read coverages for each cell on those regions to "
    "a synthetic cell with poisson distributed reads using the Jensen-Shannon "
    "distance. Cells with high enrichment of signals show a higher JSD "
    "compared to cells whose signal is homogeneously distributed."
)

USAGE = "scJSD -b sample.bam -bc barcodes.txt -o jsd.tsv"


_SAMPLING = "Sampling options"


def sampling_options() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(_SAMPLING)
    group.add_argument(
        "-n",
        "--numberOfSamples",
        dest="numberOfSamples",
        type=int,
        default=100000,
        help="The number of bins that are sampled from the genome, for which the "
        "overlapping number of reads is computed.",
    )
    group.add_argument(
        "--skipZeros",
        dest="skipZeros",
        action="store_true",
        help="If set, regions with zero overlapping reads for *all* given BAM files "
        "are ignored. This results in a reduced number of read counts compared to "
        "--numberOfSamples.",
    )
    return parser


def parse_arguments() -> argparse.ArgumentParser:
    return ParserCommon.build_parser(
        DESCRIPTION,
        USAGE,
        [
            ParserCommon.input_output_options(["bam_files", "barcodes", "out_file"]),
            ParserCommon.bam_options(
                [
                    "cell_tag",
                    "group_tag",
                    "labels",
                    "smart_labels",
                    "blacklist",
                    "chr_to_skip",
                    "bin_size",
                ],
                defaults={"bin_size": 10000},
            ),
            ParserCommon.filter_options(["duplicate_filter", "min_aligned_fraction"]),
            ParserCommon.read_options([
                "min_mapping_quality",
                "sam_flag_include",
                "sam_flag_exclude",
                "min_fragment_length",
                "max_fragment_length",
            ]),
            sampling_options(),
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
        bam_files=args.bamfiles,
        barcodes=args.barcodes,
        out_file=args.outFile,
        cell_tag=args.cellTag,
        group_tag=args.groupTag,
        labels=args.labels,
        smart_labels=args.smartLabels,
        blacklist=args.blacklist,
        chr_to_skip=args.chrToSkip,
        bin_size=args.binSize,
        duplicate_filter=args.duplicateFilter,
        min_mapping_quality=args.minMappingQuality,
        sam_flag_include=args.samFlagInclude,
        sam_flag_exclude=args.samFlagExclude,
        min_fragment_length=args.minFragmentLength,
        max_fragment_length=args.maxFragmentLength,
        min_aligned_fraction=args.minAlignedFraction,
        number_of_samples=args.numberOfSamples,
        skip_zeros=args.skipZeros,
        number_of_processors=args.numberOfProcessors,
        verbose=args.verbose,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
