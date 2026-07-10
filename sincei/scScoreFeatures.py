#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import warnings

from sincei import ParserCommon


def get_args():
    """Get scScoreFeatures-specific arguments."""
    parser = argparse.ArgumentParser(add_help=False, conflict_handler="resolve")

    scoring_opts = parser.add_argument_group("Common Options")

    scoring_opts.add_argument(
        "--features",
        "-f",
        help="Path to the BED or GTF file containing the features to use for aggregation/scoring.",
        dest="GTF",
        metavar="FILE",
        type=str,
        required=True,
    )

    scoring_opts.add_argument(
        "--overlapPolicy",
        "-op",
        help="Policy for handling regions present in .h5ad input file that only partially overlap regions present in "
        "``--features``.\n"
        "Options are: \n "
        "- ``partial``: count reads in anndata regions proportionally to the overlap fraction,"
        "(counts_considered = feature_counts * overlap_length / region_length). \n "
        "- ``all``: count all reads in the partially overlapping anndata regions. \n "
        "- ``none``: Only count reads in anndata regions that are fully contained within BED/GTF regions. \n"
        "Default: %(default)s.",
        choices=["partial", "all", "none"],
        type=str,
        default="partial",
        required=False,
    )

    scoring_opts.add_argument(
        "--centerScores",
        "-cs",
        help="If flag is set, center and scale the scores to unit variance and zero mean. Default: %(default)s.",
        action="store_true",
        required=False,
    )

    scoring_opts.add_argument(
        "--numberOfProcessors",
        "-p",
        help='Number of processors to use. Type "max/2" to '
        'use half the maximum number of processors or "max" '
        'to use all available processors. Default: "max".',
        metavar="INT",
        type=ParserCommon.numberOfProcessors,
        default=ParserCommon.numberOfProcessors("max"),
        required=False,
    )

    scoring_opts.add_argument(
        "--bedScoreFilter",
        "-bsf",
        help="Provide a range (two values separated by space), or a threshold (upper limit) of score to determine which input features to consider for scoring. "
        "Used only when the input is a BED file containing scores (stored in the 5th column). Default: %(default)s.",
        metavar="LIST",
        nargs="+",
        default=None,
        required=False,
    )

    scoring_opts.add_argument(
        "--maxRegion",
        "-mr",
        help="Maximum region size (in kb) upstream and downstream of the genes to consider for "
        "activity calculation. Default: %(default)s.",
        metavar="INT",
        type=int,
        default=100,
        required=False,
    )

    scoring_opts.add_argument(
        "--normalizeGeneLengths",
        help="Flag to indicate whether to apply length normalization to the input genes. "
        "If provided, gene scores are normalized w.r.t. gene length in the input GTF/BED file. "
        "Default: %(default)s.",
        action="store_true",
        required=False,
    )

    #    activities_opts.add_argument(
    #        "--geneBody",
    #        help="Flag to indicate whether the entire gene body is weighted as 1 (like the TSS). If provided, decay starts beyond gene body. "
    #        "By default, the weight decay starts from TSS. Default: %(default)s.",
    #        action="store_true",
    #        required=False,
    #    )

    #    activities_opts.add_argument(
    #        "--decay",
    #        "-d",
    #        help="Decay parameter for calculating distance weights. Higher values lead to faster decay "
    #        "with distance. Weights are calculated as ``exp(-decay * distance_in_kb / 10)``. "
    #        "Only used with ``--mode activities``. Default: %(default)s.",
    #        metavar="FLOAT",
    #        type=float,
    #        default=0.75,
    #        required=False,
    #    )

    #    activities_opts.add_argument(
    #        "--excludeInRange",
    #        help="Exclude regions that overlap other features from contributing to activity score of the input genes. "
    #        "This could help avoid spurious correlations between the target genes and the neighboring genes "
    #        "(in particular for promoter-enriched signals, such as H3K4me3). Options are: "
    #        "'TSS': exclude features overlapping the TSS of other genes. "
    #        "'genes': exclude features overlapping the bodies of other genes. "
    #        "Default: %(default)s.",
    #        choices=["TSS", "genes"],
    #        default=None,
    #        required=False,
    #    )

    return parser


def parseArguments(args=None):
    io_args = ParserCommon.inputOutputOptions(opts=["h5adfile", "outFile"], requiredOpts=["input", "outFile"])
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        prog="scScoreFeatures",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[io_args, get_args(), other_args],
        description="""
``scScoreFeatures`` aggregated region-level signal in a pre-processed .h5ad object into gene-level scores
based on a user-provided BED/GTF file.
""",
        usage="""
scScoreFeatures -i INPUT_binned.h5ad --features <bed or gtf> -o OUTPUT_scores.h5ad
""",
        add_help=False,
    )

    # If no arguments are provided, show help and exit
    if args is None and len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    return parser


def main(args=None):
    """Main entry point for scScoreFeatures."""

    args = parseArguments().parse_args(args)

    if not args.verbose:
        warnings.filterwarnings("ignore")
    # check scoreFilter
    if args.bedScoreFilter is not None:
        if len(args.bedScoreFilter) == 1:  # single-score: use as upper limit
            args.bedScoreFilter = [0, args.bedScoreFilter[0]]
        if len(args.bedScoreFilter) > 2:  # single-score: use as upper limit
            raise ValueError("Please provide a single value, or two values for bedScoreFilter")
        args.bedScoreFilter = [int(x) for x in args.bedScoreFilter]

    # Import here to avoid circular imports
    import anndata as ad
    from sincei.FeatureScorer import FeatureScorer

    # Load input AnnData
    adata = ParserCommon.validateAnndata(ad.read_h5ad(args.input), args.input)

    adata_out = FeatureScorer(
        adata=adata,
        gtf=args.GTF,
        mode="aggregate",
        overlap_policy=args.overlapPolicy,
        center_scores=args.centerScores,
        bedFilter=args.bedScoreFilter,
        decay=None,
        max_region=args.maxRegion,
        gene_body=True,
        gene_size_factor=args.normalizeGeneLengths,
        exclude_in_range=None,
        n_threads=args.numberOfProcessors,
        verbose=args.verbose,
    )

    if adata_out.uns is not None:
        print(".uns copied from input to output")
        adata_out.uns = adata.uns
    if adata_out.obsm is not None:
        print(".obsm copied from input to output")
        adata_out.obsm = adata.obsm
    if adata_out.obsp is not None:
        print(".obsp copied from input to output")
        adata_out.obsp = adata.obsp

    # Save output
    adata_out.write_h5ad(args.outFile)
    sys.stdout.write(f"Output saved to {args.outFile}\n")


if __name__ == "__main__":
    main()
