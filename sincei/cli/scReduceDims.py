#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import sys
import logging

from . import ParserCommon

DESCRIPTION = """
``scReduceDims`` performs dimensionality reduction on the input count matrix (output of scCountReads)
and 2D projection (UMAP) of the cells. The result is an updated h5ad object with the dimensionality
reduction results and UMAP coordinates in the ``.obsm`` field.

``scReduceDims`` provides the following dimensionality reduction methods:

* LSA: Latent Semantic Analysis.
* LDA: Latent Dirichlet Allocation.
* logPCA: Principal Component Analysis preceded by a logarithm transform.
* glmPCA: generalized PCA, with an exponential family distribution such as Poisson, Bernoulli, etc.

Each method is its own subcommand. For detailed help on each method run:

    scReduceDims <method> -h
"""

USAGE = "scReduceDims <method> -i cellCounts.h5ad -o cellCounts_pca.h5ad"

METHODS = ["LSA", "LDA", "logPCA", "glmPCA"]


def parse_arguments(args: list[str] | None = None) -> argparse.ArgumentParser:
    io_args = ParserCommon.inputOutputOptions(opts=["h5adfile", "outFile"], requiredOpts=["outFile"])
    plot_args = ParserCommon.plotOptions()
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=DESCRIPTION,
        usage=USAGE,
    )

    subparsers = parser.add_subparsers(title="Method")

    # LSA
    subparsers.add_parser(
        "LSA",
        parents=[io_args, get_args_common(), get_args_lsa(), other_args],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Use Latent Semantic Analysis (LSA) for dimensionality reduction.",
        usage="scReduceDims LSA -i cellCounts.h5ad -o cellCounts_lsa.h5ad",
        add_help=False,
    )

    # LDA
    subparsers.add_parser(
        "LDA",
        parents=[io_args, get_args_common(), get_args_lda(), other_args],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Use Latent Dirichlet Allocation (LDA) for dimensionality reduction.",
        usage="scReduceDims LDA -i cellCounts.h5ad -o cellCounts_lsa.h5ad",
        add_help=False,
    )

    # logPCA
    subparsers.add_parser(
        "logPCA",
        parents=[io_args, get_args_common(), other_args],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Use log transform + Principal Component Analysis (PCA) for dimensionality reduction.",
        usage="scReduceDims logPCA -i cellCounts.h5ad -o cellCounts_lsa.h5ad",
        add_help=False,
    )

    # logPCA
    subparsers.add_parser(
        "glmPCA",
        parents=[io_args, get_args_common(), get_args_glmpca(), other_args],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Use generalized Principal Component Analysis (glmPCA) for dimensionality reduction.",
        usage="scReduceDims glmPCA -i cellCounts.h5ad -gf poisson -o cellCounts_lsa.h5ad",
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
    parser = argparse.ArgumentParser(add_help=False)
    general = parser.add_argument_group("Clustering Options")

    general.add_argument(
        "-n",
        "--nComps",
        metavar="INT",
        default=20,
        type=int,
        help="Number of principal components or topics to reduce the dimensionality to. "
        "Use higher number for samples with more expected heterogenity. (Default: %(default)s)",
    )

    general.add_argument(
        "-nk",
        "--nNeighbors",
        metavar="INT",
        default=30,
        type=int,
        help="Number of nearest neighbours to consider for UMAP. This number should be chosen "
        "considering the total number of cells and expected number of clusters. Smaller numbers "
        "will lead to more fragmented clusters. (Default: %(default)s)",
    )

    return parser


def get_args_lsa() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    lsa_options = parser.add_argument_group("LSA Options")

    lsa_options.add_argument(
        "--binarize",
        action="store_true",
        help="Binarize the counts per region before dimensionality reduction.",
    )

    lsa_options.add_argument(
        "-om",
        "--outFileTrainedModel",
        metavar="STR",
        type=str,
        required=False,
        help="The output file for the trained LSA model. The saved model can be used later to "
        "embed/compare new cells to the existing cluster of cells.",
    )

    return parser


def get_args_lda() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    lda_options = parser.add_argument_group("LDA Options")

    lda_options.add_argument(
        "--binarize",
        action="store_true",
        help="Binarize the counts per region before dimensionality reduction (only for LSA/LDA)",
    )

    lda_options.add_argument(
        "-om",
        "--outFileTrainedModel",
        metavar="STR",
        type=str,
        required=False,
        help="The output file for the trained LDA model. The saved model can be used later to "
        "embed/compare new cells to the existing cluster of cells.",
    )

    lda_options.add_argument(
        "--nPasses",
        metavar="INT",
        type=int,
        default=5,
        help="Number of passes through the corpus for LDA model fitting. (Default: %(default)s)",
    )

    lda_options.add_argument(
        "--nIterations",
        metavar="INT",
        type=int,
        default=50,
        help="Number of iterations per pass for LDA model fitting. (Default: %(default)s)",
    )

    lda_options.add_argument(
        "--alpha",
        metavar="FLOAT",
        type=float,
        default=1.0,
        help="Prior to initialize cell-topic vectors. (Default: %(default)s)",
    )

    lda_options.add_argument(
        "--eta",
        metavar="FLOAT",
        type=float,
        default=0.1,
        help="Prior to initialize feature-topic vectors. (Default: %(default)s)",
    )

    lda_options.add_argument(
        "--gammaThreshold",
        metavar="FLOAT",
        type=float,
        nargs="?",
        help="Minimum change in the topic matrix to stop the LDA model fitting. If not given, the "
        "model is fit for the number of passes and iterations specified above.",
    )

    return parser


def get_args_glmpca() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    glmpca_options = parser.add_argument_group("glmPCA Options")

    ### PLACEHOLDER
    EXPONENTIAL_FAMILY_DICT = {
        "gaussian": "Gaussian",
        "poisson": "Poisson",
        "bernoulli": "Bernoulli",
        "beta": "Beta",
        "gamma": "Gamma",
        "lognormal": "LogNormal",
        "log_normal": "LogNormal",
        "sigmoid_beta": "SigmoidBeta",
    }

    glmpca_options.add_argument(
        "-gf",
        "--glmPCAfamily",
        type=str,
        choices=EXPONENTIAL_FAMILY_DICT.keys(),
        default="poisson",
        help="The choice of exponential family distribution to use for glmPCA method. " "(Default: %(default)s)",
    )

    return parser


def main(args: list[str] | None = None) -> int:
    args = parse_arguments().parse_args(args)

    for arg in args.__dict__:
        logging.debug(f"{arg}: {args.__dict__[arg]}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
