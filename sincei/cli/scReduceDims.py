#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys

from . import ParserCommon
from . import _parsers as backend

DESCRIPTION = (
    "Perform dimensionality reduction and UMAP projection on a cell-by-feature "
    "matrix.\n\n"
    "``scReduceDims`` performs dimensionality reduction on the input count matrix "
    "(output of scCountReads) and 2D projection (UMAP) of the cells. The result is an "
    "updated h5ad object with the dimensionality reduction results and UMAP "
    "coordinates in the ``.obsm`` field.\n\n"
    "``scReduceDims`` provides the following dimensionality reduction methods:\n"
    "* LSA: Latent Semantic Analysis.\n"
    "* LDA: Latent Dirichlet Allocation.\n"
    "* logPCA: Principal Component Analysis preceded by a logarithm transform.\n"
    "* glmPCA: generalized PCA, with an exponential family distribution such as "
    "Poisson, Bernoulli, etc."
)

USAGE = "scReduceDims {LSA,LDA,logPCA,glmPCA} -i counts.h5ad -o reduced.h5ad"

_REDUCTION = "Dimensionality reduction options"
_LDA = "LDA options"
_GLMPCA = "glmPCA options"


def reduction_options(binarize: bool) -> argparse.ArgumentParser:
    """Options shared by the dimensionality-reduction subcommands.

    ``binarize`` selects the two options that only LSA and LDA carry.
    """
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(_REDUCTION)

    group.add_argument(
        "-n",
        "--nComps",
        dest="nComps",
        type=int,
        default=20,
        help="Number of principal components or topics to reduce the dimensionality "
        "to. Use a higher number for samples with more expected heterogeneity.",
    )
    group.add_argument(
        "-nk",
        "--nNeighbors",
        dest="nNeighbors",
        type=int,
        default=30,
        help="Number of nearest neighbours to consider for UMAP. Choose this "
        "considering the total number of cells and the expected number of clusters; "
        "smaller numbers lead to more fragmented clusters.",
    )
    if binarize:
        group.add_argument(
            "--binarize",
            action="store_true",
            help="Binarize the counts per region before dimensionality reduction.",
        )
        group.add_argument(
            "-om",
            "--outFileTrainedModel",
            dest="outFileTrainedModel",
            default=None,
            help="The output file for the trained model. The saved model can be used "
            "later to embed/compare new cells to the existing cluster of cells.",
        )
    return parser


def lda_options() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(_LDA)

    group.add_argument(
        "--nPasses",
        dest="nPasses",
        type=int,
        default=5,
        help="Number of passes through the corpus for LDA model fitting.",
    )
    group.add_argument(
        "--nIterations",
        dest="nIterations",
        type=int,
        default=50,
        help="Number of iterations per pass for LDA model fitting.",
    )
    group.add_argument(
        "--alpha",
        type=float,
        default=1.0,
        help="Prior to initialize cell-topic vectors.",
    )
    group.add_argument(
        "--eta",
        type=float,
        default=0.1,
        help="Prior to initialize feature-topic vectors.",
    )
    group.add_argument(
        "--gammaThreshold",
        dest="gammaThreshold",
        type=float,
        default=None,
        help="Minimum change in the topic matrix to stop the LDA model fitting. If not "
        "given, the model is fit for the number of passes and iterations specified "
        "above.",
    )
    return parser


def glmpca_options() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(_GLMPCA)

    group.add_argument(
        "-gf",
        "--glmPCAfamily",
        dest="glmPCAfamily",
        metavar="FAMILY",
        choices=ParserCommon.GLMPCA_FAMILIES,
        default="poisson",
        help="The choice of exponential family distribution to use for the glmPCA "
        "method.",
    )
    return parser


_SUBCOMMANDS = {
    "LSA": ("Use Latent Semantic Analysis (LSA).", True, []),
    "LDA": ("Use Latent Dirichlet Allocation (LDA).", True, [lda_options]),
    "logPCA": ("Use log transform + PCA.", False, []),
    "glmPCA": ("Use generalized PCA.", False, [glmpca_options]),
}


def parse_arguments() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=DESCRIPTION,
        usage=USAGE,
        add_help=False,
    )
    parser.add_argument(
        "-h", "--help", action="help", help="Show this message and exit."
    )
    subparsers = parser.add_subparsers(dest="method", metavar="{LSA,LDA,logPCA,glmPCA}")

    for name, (help_text, binarize, extra) in _SUBCOMMANDS.items():
        sub = subparsers.add_parser(
            name,
            parents=[
                ParserCommon.input_output_options(["h5ad_file", "out_file"]),
                reduction_options(binarize),
                *(factory() for factory in extra),
                ParserCommon.other_options(),
            ],
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=help_text,
            help=help_text,
            add_help=False,
            conflict_handler="resolve",
        )
        sub.set_defaults(method=name)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = parse_arguments()
    if argv is None and len(sys.argv) == 1:
        parser.print_help()
        return 0
    args = parser.parse_args(argv)
    if not getattr(args, "method", None):
        parser.print_help()
        return 0

    backend.log_parameters(
        input=args.input,
        out_file=args.outFile,
        n_comps=args.nComps,
        n_neighbors=args.nNeighbors,
        number_of_processors=args.numberOfProcessors,
        verbose=args.verbose,
    )
    if args.method in ("LSA", "LDA"):
        backend.log_parameters(
            binarize=args.binarize,
            out_file_trained_model=args.outFileTrainedModel,
        )
    if args.method == "LDA":
        backend.log_parameters(
            n_passes=args.nPasses,
            n_iterations=args.nIterations,
            alpha=args.alpha,
            eta=args.eta,
            gamma_threshold=args.gammaThreshold,
        )
    if args.method == "glmPCA":
        backend.log_parameters(glmpca_family=args.glmPCAfamily)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
