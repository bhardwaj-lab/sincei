#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
import warnings

from sincei import ParserCommon

logger = logging.getLogger()


def get_args():
    """Get scScoreFeatures-specific arguments."""
    parser = argparse.ArgumentParser(add_help=False, conflict_handler="resolve")

    scoring_opts = parser.add_argument_group("Scoring Options")

    scoring_opts.add_argument(
        "--mode",
        "-m",
        help="When in 'activities' mode, calculates weighted gene activity scores using exponential "
        "decay from gene body/TSS. \n"
        "'VCR' mode calculates simple sum of counts within gene body. "
        "Only used with --GTF. Required when using --GTF.",
        choices=["activities", "VCR"],
        default=None,
        required=False,
    )

    scoring_opts.add_argument(
        "--VCR",
        "-VCR",
        help="Path to the BED file containing the variable chromatin regions (VCRs) to use for count aggregation.",
        metavar="VCRs.BED",
        type=str,
        required=False,
    )

    scoring_opts.add_argument(
        "--penalty",
        "-pen",
        help="Penalty value in the VCR BED file (5th column) to determine which VCRs to use for aggregation. "
        "Only used with --VCR. Default: %(default)s.",
        metavar="FLOAT",
        type=float,
        default=0.05,
        required=False,
    )

    scoring_opts.add_argument(
        "--GTF",
        "-GTF",
        help="Path to the GTF file containing gene annotations for which to compute gene activity scores.",
        metavar="GTF",
        type=str,
        required=False,
    )

    scoring_opts.add_argument(
        "--decay",
        "-d",
        help="Decay parameter for calculating distance weights. Higher values lead to faster decay "
        " with distance. Weights are calculated as `exp(-decay * distance_in_kb)`. "
        "Only used with --GTF in activities mode. Default: %(default)s.",
        metavar="FLOAT",
        type=float,
        default=0.75,
        required=False,
    )

    scoring_opts.add_argument(
        "--maxRegion",
        "-mr",
        help="Maximum region size (bp) around the gene (upstream and downstream) to consider for "
        "gene activity calculation. Only used with --GTF. Default: %(default)s.",
        metavar="INT",
        type=int,
        default=200000,
        required=False,
    )

    scoring_opts.add_argument(
        "--geneBody",
        help="Whether the gene body is weighted as 1 (like TSS). If True, decay starts beyond gene body. "
        "If False, decay starts from TSS. Only used with --GTF. Default: %(default)s.",
        action="store_true",
        default=True,
        required=False,
    )

    scoring_opts.add_argument(
        "--geneSizeFactor",
        help="Apply gene length normalization factor. Only used with --GTF. Default: %(default)s.",
        action="store_true",
        default=True,
        required=False,
    )

    scoring_opts.add_argument(
        "--excludeInRange",
        help="Exclude regions of other genes from contributing to gene activity scores. "
        "'TSS': exclude TSS ± extendTSS regions. 'genes': exclude gene bodies extended upstream by extendTSS. "
        "None: no exclusion. Only used with --GTF. Default: %(default)s.",
        choices=["TSS", "genes"],
        default=None,
        required=False,
    )

    scoring_opts.add_argument(
        "--extendTSS",
        help="Number of base pairs to extend around TSS for exclusion regions when using --excludeInRange. "
        "Only used with --GTF. Default: %(default)d.",
        metavar="INT",
        type=int,
        default=20,
        required=False,
    )

    common_opts = parser.add_argument_group("Common Options")

    common_opts.add_argument(
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

    common_opts.add_argument(
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
    
    # Score gene activities from chromatin features
    scScoreFeatures -i chrom_features.h5ad --GTF genome.gtf -o gene_activities.h5ad
    
    # Aggregate chromatin bins into VCRs
    scScoreFeatures -i chrom_bins.h5ad --VCR VCRs.bed --penalty 0.05 -o chrom_VCRs.h5ad
""",
        usage="""
scScoreFeatures -m activities -i INPUT.h5ad --GTF genome.gtf -o OUTPUT.h5ad [options]
scScoreFeatures -m VCR -i INPUT.h5ad --VCR VCRs.bed -o OUTPUT.h5ad [options]
""",
        add_help=False,
    )

    # If no arguments are provided, show help and exit
    if args is None and len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args(args)

    # Validate that either --GTF or --VCR is provided (but not both)
    if args.GTF and args.VCR:
        parser.error("Cannot specify both --GTF and --VCR.")
    if not args.GTF and not args.VCR:
        parser.error("Must specify either --GTF or --VCR.")

    # Validate mode is specified
    if args.mode is None:
        parser.error("You must specify --mode. Choose either 'activities' or 'VCR'.")

    # Warn --penalty is only used with --VCR
    if args.penalty != 0.05 and not args.VCR and args.verbose:
        sys.stderr.write("Ignoring --penalty. It can only be used with --VCR (aggregation mode).\n")

    return args


def main(args=None):
    """Main entry point for scScoreFeatures."""

    args = parse_arguments(args)

    if not args.verbose:
        logger.setLevel(logging.CRITICAL)
        warnings.filterwarnings("ignore")

    # Import here to avoid circular imports
    import anndata as ad
    from sincei.feature_scorer import feature_scorer

    # Load input AnnData
    adata = ad.read_h5ad(args.input)

    # Determine mode based on which arguments are present
    if args.VCR:
        # VCR aggregation mode - use VCR BED file
        adata_out = feature_scorer(
            adata=adata,
            gtf=args.VCR,
            mode="VCR",
            penalty=args.penalty,
            decay=args.decay,
            max_region=args.maxRegion,
            gene_body=args.geneBody,
            gene_size_factor=args.geneSizeFactor,
            exclude_in_range=args.excludeInRange,
            extend_TSS=args.extendTSS,
            chrs_to_skip=args.chrToSkip,
            n_threads=args.numberOfProcessors,
        )
    elif args.GTF:
        # Gene scoring mode - use GTF with user-specified mode and parameters
        adata_out = feature_scorer(
            adata=adata,
            gtf=args.GTF,
            mode=args.mode,
            decay=args.decay,
            max_region=args.maxRegion,
            gene_body=args.geneBody,
            gene_size_factor=args.geneSizeFactor,
            exclude_in_range=args.excludeInRange,
            extend_TSS=args.extendTSS,
            chrs_to_skip=args.chrToSkip,
            n_threads=args.numberOfProcessors,
        )

    # Save output
    adata_out.write_h5ad(args.outFile)
    sys.stdout.write(f"Output saved to {args.outFile}\n")


if __name__ == "__main__":
    main()
