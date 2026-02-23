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
        help="When in ``activities`` mode, calculates weighted gene activity scores using exponential "
        "decay per cell for each gene body/TSS or region. \n"
        "``aggregate`` mode calculates the simple sum of counts per cell for each gene body or region.",
        choices=["aggregate", "activities"],
        default=None,
        required=True,
    )

    scoring_opts.add_argument(
        "--GTF",
        "-GTF",
        "--BED",
        "-BED",
        help="Path to the BED/GTF file containing the regions to use for aggregation/feature scoring.",
        dest="GTF",
        metavar="FILE",
        type=str,
        required=True,
    )

    scoring_opts.add_argument(
        "--overlapPolicy",
        "-op",
        help="Policy for handling adata features that only partially overlap regions in the BED/GTF provided.\n"
        " Options are:\n"
        "    - ``partial``: count reads in anndata feature proportionally to the overlap fraction, \n"
        "       read as counts_considered = feature_counts * overlap_length / region_length.\n"
        "    - ``all``: count all reads in the partially overlapping anndata feature.\n"
        "    - ``none``: exclude reads from partially overlapping anndata features, in other words, only\n"
        "      count reads in anndata features fully contained within BED/GTF regions.\n"
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
        help="Penalty value to determine which VCRs to use for aggregation for a VCR BED file with "
        "multiple penalties (5th column). Default: %(default)s.",
        metavar="FLOAT",
        type=float,
        default=None,
        required=False,
    )

    activities_opts.add_argument(
        "--decay",
        "-d",
        help="Decay parameter for calculating distance weights. Higher values lead to faster decay "
        " with distance. Weights are calculated as ``exp(-decay * distance_in_kb)``. "
        "Only used with --GTF in activities mode. Default: %(default)s.",
        metavar="FLOAT",
        type=float,
        default=0.75,
        required=False,
    )

    activities_opts.add_argument(
        "--maxRegion",
        "-mr",
        help="Maximum region size (bp) around the gene (upstream and downstream) to consider for "
        "gene activity calculation. Default: %(default)s.",
        metavar="INT",
        type=int,
        default=200000,
        required=False,
    )

    activities_opts.add_argument(
        "--geneBody",
        help="Whether the gene body is weighted as 1 (like TSS). If True, decay starts beyond gene body. "
        "If False, decay starts from TSS. Default: %(default)s.",
        action="store_true",
        default=True,
        required=False,
    )

    activities_opts.add_argument(
        "--geneSizeFactor",
        help="Apply gene length normalization factor. Default: %(default)s.",
        action="store_true",
        default=True,
        required=False,
    )

    activities_opts.add_argument(
        "--excludeInRange",
        help="Exclude regions of other genes from contributing to gene activity scores. "
        "'TSS': exclude TSS ± extendTSS regions. 'genes': exclude gene bodies extended upstream by extendTSS. "
        "None: no exclusion. Default: %(default)s.",
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
    from sincei.feature_scorer import feature_scorer

    # Load input AnnData
    adata = ad.read_h5ad(args.input)

    adata_out = feature_scorer(
        adata=adata,
        gtf=args.GTF,
        mode=args.mode,
        overlap_policy=args.overlapPolicy,
        penalty=args.penalty,
        decay=args.decay,
        max_region=args.maxRegion,
        gene_body=args.geneBody,
        gene_size_factor=args.geneSizeFactor,
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
