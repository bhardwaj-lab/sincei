#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys

from sincei import _sincei as internal

from . import ParserCommon
from . import _parsers as backend

DESCRIPTION = (
    "Aggregate a binned chromatin count matrix into per-feature scores.\n\n"
    "``scScoreFeatures`` sums the counts of the bins overlapping each feature in "
    "``--features``, producing a cells x features matrix. The features can be genes "
    "(from a GTF/GFF) or Variable Chromatin Regions (from a BED file produced by "
    "scFindVCRs)."
)

USAGE = "scScoreFeatures -i counts.h5ad --features genes.gtf -o scores.h5ad"

_SCORING = "Scoring options"


def scoring_options() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(_SCORING)

    group.add_argument(
        "--features",
        metavar=".bed/.gtf/.gff",
        required=True,
        help="Path to the BED, GTF or GFF file containing the features to score.",
    )
    group.add_argument(
        "-op",
        "--overlapPolicy",
        dest="overlapPolicy",
        metavar="POLICY",
        choices=ParserCommon.OVERLAP_POLICIES,
        default="partial",
        help="How to treat a bin in the .h5ad input that only partially overlaps a "
        "region in --features.\n\n"
        "partial: count the fraction of the bin lying inside the region; "
        "all: count the whole bin; "
        "none: ignore the bin unless it lies wholly inside the region.",
    )
    group.add_argument(
        "-pen",
        "--penalty",
        dest="penalty",
        type=float,
        default=None,
        help="Penalty value to determine which VCRs to score. Used only when the input "
        "is a BED file created with ``scFindVCRs`` with a range of penalties (stored "
        "in the 5th column).",
    )
    return parser


def parse_arguments() -> argparse.ArgumentParser:
    return ParserCommon.build_parser(
        DESCRIPTION,
        USAGE,
        [
            ParserCommon.input_output_options(["h5ad_file", "out_file"]),
            scoring_options(),
            ParserCommon.gtf_gff_options(
                ["transcript_id", "exon_id", "transcript_id_tag", "metagene"],
                metagene_flags=("-m", "--metagene"),
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

    if args.verbose:
        backend.log_parameters(
            input=args.input,
            out_file=args.outFile,
            features=args.features,
            overlap_policy=args.overlapPolicy,
            penalty=args.penalty,
            feature_type=args.transcriptID,
            exon_id=args.exonID,
            name_attr=args.transcriptIDtag,
            metagene=args.metagene,
            number_of_processors=args.numberOfProcessors,
        )
    result_path = internal.score_features(
        args.input,
        args.features,
        args.outFile,
        overlap_policy=args.overlapPolicy,
        penalty=args.penalty,
        feature_type=args.transcriptID,
        exon_type=args.exonID,
        name_attr=args.transcriptIDtag,
        metagene=args.metagene,
        num_threads=args.numberOfProcessors,
    )
    print(result_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
