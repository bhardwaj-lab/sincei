#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
import warnings

from sincei import ParserCommon


def get_args():
    """Get scScoreFeatures-specific arguments."""
    parser = argparse.ArgumentParser(add_help=False, conflict_handler="resolve")

    scoring_opts = parser.add_argument_group("Common Options")
    aggregate_opts = parser.add_argument_group("Aggregate Mode Options")
    activities_opts = parser.add_argument_group("Activities Mode Options")

    scoring_opts.add_argument(
        "--mode",
        "-m",
        help="The ``activities`` mode calculates weighted gene activity scores per cell using exponential "
        "decay of signal around each input feature (such as gene-body, TSS, or region). \n"
        "``aggregate`` mode calculates a simple sum of counts per cell for each input feature.",
        choices=["aggregate", "activities"],
        default=None,
        required=True,
    )

    scoring_opts.add_argument(
        "--features",
        help="Path to the .BED or .GTF formatted file containing the features to use for aggregation/scoring.",
        dest="GTF",
        metavar="FILE",
        type=str,
        required=True,
    )

    scoring_opts.add_argument(
        "--overlapPolicy",
        "-op",
        help="Policy for handling regions present in .h5ad input file that only partially overlap regions present in --features.\n"
        " Options are:\n"
        "    - ``partial``: count reads in anndata regions proportionally to the overlap fraction, \n"
        "       read as counts_considered = feature_counts * overlap_length / region_length.\n"
        "    - ``all``: count all reads in the partially overlapping anndata regions.\n"
        "    - ``none``: exclude reads from partially overlapping anndata regions, in other words, only\n"
        "      count reads in anndata regions fully contained within BED/GTF regions.\n"
        "Default is %(default)s.",
        choices=["partial", "all", "none"],
        type=str,
        default="partial",
        required=False,
    )

    scoring_opts.add_argument(
        "--chrToSkip",
        help="List of chromosome names to skip from the analysis. "
        "Regions on these chromosomes will be excluded. "
        "Useful for skipping mitochondrial, X chromosome, or unplaced contigs. "
        "Multiple chromosomes can be specified, e.g. ``--chrToSkip chrM chrX``.",
        metavar="CHR",
        nargs="+",
        default=None,
        required=False,
    )

    scoring_opts.add_argument(
        "--numberOfProcessors",
        "-p",
        help='Number of processors to use. Type "max/2" to '
        'use half the maximum number of processors or "max" '
        'to use all available processors. Default: "max"',
        metavar="INT",
        type=ParserCommon.numberOfProcessors,
        default=ParserCommon.numberOfProcessors("max"),
        required=False,
    )

    aggregate_opts.add_argument(
        "--penalty",
        "-pen",
        help="Penalty value to determine which VCRs to use for aggregation. Used only when the input is "
        "a BED file created with `scFindVCRs` with a range of penalties (stored in the 5th column). Default: %(default)s.",
        metavar="FLOAT",
        type=float,
        default=None,
        required=False,
    )

    activities_opts.add_argument(
        "--decay",
        "-d",
        help="Decay parameter for calculating distance weights. Higher values lead to faster decay "
        "with distance. Weights are calculated as ``exp(-decay * distance_in_kb)``. "
        "Only used with `--mode activities` . Default: %(default)s.",
        metavar="FLOAT",
        type=float,
        default=0.75,
        required=False,
    )

    activities_opts.add_argument(
        "--maxRegion",
        "-mr",
        help="Maximum region size (in kb) upstream and downstream of the genes to consider for "
        "activity calculation. Default: %(default)s.",
        metavar="INT",
        type=int,
        default=200,
        required=False,
    )

    activities_opts.add_argument(
        "--geneBody",
        help="Flag to indicate whether the entire gene body is weighted as 1 (like TSS). If provided, decay starts beyond gene body. "
        "By default, the weight decay starts from TSS. Default: %(default)s.",
        action="store_true",
        required=False,
    )

    activities_opts.add_argument(
        "--normalizeGeneLengths",
        help="Flag to indicate whether to apply length normalization to the input genes. "
        "If provided. Each gene is normalized w.r.t the average gene length in the input GTF/BED file "
        "Default: %(default)s.",
        action="store_true",
        required=False,
    )

    activities_opts.add_argument(
        "--excludeInRange",
        help="Exclude regions that overlap other features from contributing to activity score of the input genes. "
        "This could help avoid spurious correlations between the target genes and the neighboring genes (in perticular for promoter-enriched ) "
        "signals, such as H3K4me3. Options are: "
        "'TSS': exclude TSS ± extendTSS regions. "
        "'genes': exclude gene bodies extended upstream by extendTSS. "
        "Default: %(default)s.",
        choices=["TSS", "genes"],
        default=None,
        required=False,
    )

    activities_opts.add_argument(
        "--extendTSS",
        help="Number of base pairs to extend around TSS for exclusion regions when using --excludeInRange. "
        "Default: %(default)d.",
        metavar="INT",
        type=int,
        default=20,
        required=False,
    )

    return parser


def parse_arguments(args=None):
    io_args = ParserCommon.inputOutputOptions(opts=["h5adfile", "outFile"], requiredOpts=["input", "outFile"])
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        prog="scScoreFeatures",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[io_args, get_args(), other_args],
        description="""
``scScoreFeatures`` computes gene activity scores from chromatin data (use --GTF) or
aggregates binned chromatin data into Variable Chromatin Regions (use --VCR with output from scFindVCRs).

Examples:
    # Aggregate chromatin bins into VCRs
    scScoreFeatures -m aggregate -i chrom_bins.h5ad --VCR VCRs.bed --penalty 0.05 -o chrom_VCRs.h5ad

    # Score gene activities from chromatin features
    scScoreFeatures -m activities -i chrom_features.h5ad --GTF genome.gtf -o gene_activities.h5ad
""",
        usage="""
scScoreFeatures -m aggregate -i INPUT.h5ad --GTF VCRs.bed -o OUTPUT.h5ad [options]
scScoreFeatures -m activities -i INPUT.h5ad --GTF genome.gtf -o OUTPUT.h5ad [options]
""",
        add_help=False,
    )

    # If no arguments are provided, show help and exit
    if args is None and len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args(args)

    # Validate mode is specified
    if args.mode is None:
        parser.error("You must specify --mode. Choose either 'activities' or 'aggregate'.")

    return args


def main(args=None):
    """Main entry point for scScoreFeatures."""

    args = parse_arguments(args)

    if not args.verbose:
        warnings.filterwarnings("ignore")

    # Import here to avoid circular imports
    import anndata as ad
    from sincei.FeatureScorer import FeatureScorer

    # Load input AnnData
    adata = ad.read_h5ad(args.input)

    adata_out = FeatureScorer(
        adata=adata,
        gtf=args.GTF,
        mode=args.mode,
        overlap_policy=args.overlapPolicy,
        penalty=args.penalty,
        decay=args.decay,
        max_region=args.maxRegion,
        gene_body=args.geneBody,
        gene_size_factor=args.normalizeGeneLengths,
        exclude_in_range=args.excludeInRange,
        extend_TSS=args.extendTSS,
        chrs_to_skip=args.chrToSkip,
        n_threads=args.numberOfProcessors,
        verbose=args.verbose,
    )

    # Save output
    adata_out.write_h5ad(args.outFile)
    sys.stdout.write(f"Output saved to {args.outFile}\n")


if __name__ == "__main__":
    main()
