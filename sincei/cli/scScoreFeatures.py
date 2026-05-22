#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import sys
import logging

from . import ParserCommon

# from sincei import _sincei

DESCRIPTION = """
``scScoreFeatures`` computes gene activity scores from chromatin data or aggregates
binned chromatin data into Variable Chromatin Regions (use output from scFindVCRs).

Detailed help for each sub-command is available by typing:

    scScoreFeatures activities -h
    scScoreFeatures aggregate -h
"""

USAGE = "scScoreFeatures <method> -i binned_data.h5ad --features features.bed  -o counts.h5ad"

METHODS = ["activities", "aggregate"]


def parse_arguments(args: list[str] | None = None) -> argparse.ArgumentParser:
    io_args = ParserCommon.inputOutputOptions(opts=["h5adfile", "outFile"], requiredOpts=["input", "outFile"])
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        prog="scScoreFeatures",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[io_args, get_args_common(), other_args],
        description=DESCRIPTION,
        usage=USAGE,
        add_help=False,
    )

    subparsers = parser.add_subparsers(title="Method")

    # Sum
    subparsers.add_parser(
        "activities",
        parents=[io_args, get_args_common(), get_args_activities(), other_args],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Calculates gene activity scores from binned chromatin data.",
        usage="scScoreFeatures -m activities -i INPUT_binned.h5ad --features genes.gtf -o OUTPUT_activities.h5ad",
        add_help=False,
    )

    # Aggregate
    subparsers.add_parser(
        "aggregate",
        parents=[io_args, get_args_common(), get_args_aggregate(), other_args],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Aggregate binned chromatin data into VCRs using the output of scFindVCRs.",
        usage="scScoreFeatures aggregate -i INPUT_binned.h5ad --features VCRs.bed  -o OUTPUT_aggregate.h5ad",
        add_help=False,
    )

    # If no arguments are provided, show help and exit
    if args is None and len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    # If only subcommand is provided, show subcommand help
    if sys.argv[1] in METHODS and len(sys.argv) == 2:
        parser.parse_args([sys.argv[1], "-h"])

    return parser


def get_args_common() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False, conflict_handler="resolve")

    general = parser.add_argument_group("Common Options")

    general.add_argument(
        "--features",
        metavar=[".bed", ".gtf"],
        help="Path to the BED or GTF file containing the features to use for aggregation/scoring.",
        dest="GTF",
        type=str,
        required=True,
    )

    general.add_argument(
        "--overlapPolicy",
        "-op",
        help="Policy for handling regions present in .h5ad input file that only partially overlap regions present in --features.\n"
        " Options are: \n "
        "    - ``partial``: count reads in anndata regions proportionally to the overlap fraction, \n "
        "       read as counts_considered = feature_counts * overlap_length / region_length. \n "
        "    - ``all``: count all reads in the partially overlapping anndata regions.\n"
        "    - ``none``: exclude reads from partially overlapping anndata regions, in other words, only \n "
        "      count reads in anndata regions fully contained within BED/GTF regions.\n"
        "Default: %(default)s.",
        choices=["partial", "all", "none"],
        type=str,
        default="partial",
        required=False,
    )

    general.add_argument(
        "--centerScores",
        "-cs",
        help="If flag is set, center and scale the scores to unit variance and zero mean. Default: %(default)s.",
        action="store_true",
        required=False,
    )

    return parser


def get_args_activities() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False, conflict_handler="resolve")
    activities_options = parser.add_argument_group("Activities Options")

    activities_options.add_argument(
        "-d",
        "--decay",
        help="Decay parameter for calculating distance weights. Higher values lead to faster decay "
        "with distance. Weights are calculated as ``exp(-decay * distance_in_kb / 10)``. "
        "Only used with ``--mode activities``. Default: %(default)s.",
        metavar="FLOAT",
        type=float,
        default=0.75,
        required=False,
    )

    activities_options.add_argument(
        "-mr",
        "--maxRegion",
        help="Maximum region size (in kb) upstream and downstream of the genes to consider for "
        "activity calculation. Default: %(default)s.",
        metavar="INT",
        type=int,
        default=100,
        required=False,
    )

    activities_options.add_argument(
        "--geneBody",
        help="Flag to indicate whether the entire gene body is weighted as 1 (like the TSS). If provided, decay starts beyond gene body. "
        "By default, the weight decay starts from TSS. Default: %(default)s.",
        action="store_true",
        required=False,
    )

    activities_options.add_argument(
        "--normalizeGeneLengths",
        help="Flag to indicate whether to apply length normalization to the input genes. "
        "If provided, gene scores are normalized w.r.t. gene length in the input GTF/BED file. "
        "Default: %(default)s.",
        action="store_true",
        required=False,
    )

    activities_options.add_argument(
        "--excludeInRange",
        help="Exclude regions that overlap other features from contributing to activity score of the input genes. "
        "This could help avoid spurious correlations between the target genes and the neighboring genes "
        "(in particular for promoter-enriched signals, such as H3K4me3). Options are: "
        "'TSS': exclude features overlapping the TSS of other genes. "
        "'genes': exclude features overlapping the bodies of other genes. "
        "Default: %(default)s.",
        choices=["TSS", "genes"],
        default=None,
        required=False,
    )

    return parser


def get_args_aggregate() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False, conflict_handler="resolve")
    aggregate_options = parser.add_argument_group("Aggregate Options")

    aggregate_options.add_argument(
        "-pen",
        "--penalty",
        help="Penalty value to determine which VCRs to use for aggregation. Used only when the input is "
        "a BED file created with ``scFindVCRs`` with a range of penalties (stored in the 5th column). Default: %(default)s.",
        metavar="FLOAT",
        type=float,
        default=None,
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
