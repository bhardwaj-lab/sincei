from __future__ import annotations

import typer

from .shared import log_parameters, normalize_processors, version_option

DESCRIPTION = """
Perform dimensionality reduction and UMAP projection on a cell-by-feature matrix.
"""


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help=DESCRIPTION,
    context_settings={"help_option_names": ["-h", "--help"]},
)

lsa_app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Use Latent Semantic Analysis (LSA).",
    context_settings={"help_option_names": ["-h", "--help"]},
)
lda_app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Use Latent Dirichlet Allocation (LDA).",
    context_settings={"help_option_names": ["-h", "--help"]},
)
logpca_app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Use log transform + PCA.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
glmpca_app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Use generalized PCA.",
    context_settings={"help_option_names": ["-h", "--help"]},
)

app.add_typer(lsa_app, name="LSA")
app.add_typer(lda_app, name="LDA")
app.add_typer(logpca_app, name="logPCA")
app.add_typer(glmpca_app, name="glmPCA")


@app.callback(invoke_without_command=True)
def main(ctx: typer.Context) -> int:
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
        raise typer.Exit()
    return 0


def _common(
    input: str,
    outFile: str,
    nComps: int,
    nNeighbors: int,
    numberOfProcessors: str,
    verbose: bool,
) -> None:
    log_parameters(
        input=input,
        outFile=outFile,
        nComps=nComps,
        nNeighbors=nNeighbors,
        numberOfProcessors=numberOfProcessors,
        verbose=verbose,
    )


@version_option("scReduceDims LSA")
@lsa_app.callback(invoke_without_command=True)
def lsa(
    input: str = typer.Option(
        ..., "-i", "--input", metavar=".h5ad", help="sincei-generated input file in .h5ad format."
    ),
    outFile: str = typer.Option(..., "-o", "--outFile", help="The file to write results to."),
    nComps: int = typer.Option(
        20,
        "-n",
        "--nComps",
        help="Number of principal components or topics to reduce the dimensionality to. Use higher number for samples with more expected heterogenity. (Default: %(default)s)",
    ),
    nNeighbors: int = typer.Option(
        30,
        "-nk",
        "--nNeighbors",
        help="Number of nearest neighbours to consider for UMAP. This number should be chosen considering the total number of cells and expected number of clusters. Smaller numbers will lead to more fragmented clusters. (Default: %(default)s)",
    ),
    binarize: bool = typer.Option(
        False, "--binarize", help="Binarize the counts per region before dimensionality reduction."
    ),
    outFileTrainedModel: str = typer.Option(
        None,
        "-om",
        "--outFileTrainedModel",
        help="The output file for the trained LSA model. The saved model can be used later to embed/compare new cells to the existing cluster of cells.",
    ),
    numberOfProcessors: str = typer.Option(
        "max",
        "-p",
        "--numberOfProcessors",
        callback=normalize_processors,
        help='Number of processors to use. Type "max/2" to use half the maximum number of processors or "max" to use all available processors. (Default: "max")',
    ),
    verbose: bool = typer.Option(False, "-v", "--verbose", help="Set to see processing messages."),
) -> int:
    _common(input, outFile, nComps, nNeighbors, numberOfProcessors, verbose)
    log_parameters(binarize=binarize, outFileTrainedModel=outFileTrainedModel)
    return 0


@version_option("scReduceDims LDA")
@lda_app.callback(invoke_without_command=True)
def lda(
    input: str = typer.Option(
        ..., "-i", "--input", metavar=".h5ad", help="sincei-generated input file in .h5ad format."
    ),
    outFile: str = typer.Option(..., "-o", "--outFile", help="The file to write results to."),
    nComps: int = typer.Option(
        20,
        "-n",
        "--nComps",
        help="Number of principal components or topics to reduce the dimensionality to. Use higher number for samples with more expected heterogenity. (Default: %(default)s)",
    ),
    nNeighbors: int = typer.Option(
        30,
        "-nk",
        "--nNeighbors",
        help="Number of nearest neighbours to consider for UMAP. This number should be chosen considering the total number of cells and expected number of clusters. Smaller numbers will lead to more fragmented clusters. (Default: %(default)s)",
    ),
    binarize: bool = typer.Option(
        False, "--binarize", help="Binarize the counts per region before dimensionality reduction (only for LSA/LDA)"
    ),
    outFileTrainedModel: str = typer.Option(
        None,
        "-om",
        "--outFileTrainedModel",
        help="The output file for the trained LDA model. The saved model can be used later to embed/compare new cells to the existing cluster of cells.",
    ),
    nPasses: int = typer.Option(
        5, "--nPasses", help="Number of passes through the corpus for LDA model fitting. (Default: %(default)s)"
    ),
    nIterations: int = typer.Option(
        50, "--nIterations", help="Number of iterations per pass for LDA model fitting. (Default: %(default)s)"
    ),
    alpha: float = typer.Option(1.0, "--alpha", help="Prior to initialize cell-topic vectors. (Default: %(default)s)"),
    eta: float = typer.Option(0.1, "--eta", help="Prior to initialize feature-topic vectors. (Default: %(default)s)"),
    gammaThreshold: float = typer.Option(
        None,
        "--gammaThreshold",
        help="Minimum change in the topic matrix to stop the LDA model fitting. If not given, the model is fit for the number of passes and iterations specified above.",
    ),
    numberOfProcessors: str = typer.Option(
        "max",
        "-p",
        "--numberOfProcessors",
        callback=normalize_processors,
        help='Number of processors to use. Type "max/2" to use half the maximum number of processors or "max" to use all available processors. (Default: "max")',
    ),
    verbose: bool = typer.Option(False, "-v", "--verbose", help="Set to see processing messages."),
) -> int:
    _common(input, outFile, nComps, nNeighbors, numberOfProcessors, verbose)
    log_parameters(
        binarize=binarize,
        outFileTrainedModel=outFileTrainedModel,
        nPasses=nPasses,
        nIterations=nIterations,
        alpha=alpha,
        eta=eta,
        gammaThreshold=gammaThreshold,
    )
    return 0


@version_option("scReduceDims logPCA")
@logpca_app.callback(invoke_without_command=True)
def logpca(
    input: str = typer.Option(
        ..., "-i", "--input", metavar=".h5ad", help="sincei-generated input file in .h5ad format."
    ),
    outFile: str = typer.Option(..., "-o", "--outFile", help="The file to write results to."),
    nComps: int = typer.Option(
        20,
        "-n",
        "--nComps",
        help="Number of principal components or topics to reduce the dimensionality to. Use higher number for samples with more expected heterogenity. (Default: %(default)s)",
    ),
    nNeighbors: int = typer.Option(
        30,
        "-nk",
        "--nNeighbors",
        help="Number of nearest neighbours to consider for UMAP. This number should be chosen considering the total number of cells and expected number of clusters. Smaller numbers will lead to more fragmented clusters. (Default: %(default)s)",
    ),
    numberOfProcessors: str = typer.Option(
        "max",
        "-p",
        "--numberOfProcessors",
        callback=normalize_processors,
        help='Number of processors to use. Type "max/2" to use half the maximum number of processors or "max" to use all available processors. (Default: "max")',
    ),
    verbose: bool = typer.Option(False, "-v", "--verbose", help="Set to see processing messages."),
) -> int:
    _common(input, outFile, nComps, nNeighbors, numberOfProcessors, verbose)
    return 0


@version_option("scReduceDims glmPCA")
@glmpca_app.callback(invoke_without_command=True)
def glmpca(
    input: str = typer.Option(
        ..., "-i", "--input", metavar=".h5ad", help="sincei-generated input file in .h5ad format."
    ),
    outFile: str = typer.Option(..., "-o", "--outFile", help="The file to write results to."),
    nComps: int = typer.Option(
        20,
        "-n",
        "--nComps",
        help="Number of principal components or topics to reduce the dimensionality to. Use higher number for samples with more expected heterogenity. (Default: %(default)s)",
    ),
    nNeighbors: int = typer.Option(
        30,
        "-nk",
        "--nNeighbors",
        help="Number of nearest neighbours to consider for UMAP. This number should be chosen considering the total number of cells and expected number of clusters. Smaller numbers will lead to more fragmented clusters. (Default: %(default)s)",
    ),
    glmPCAfamily: str = typer.Option(
        "poisson",
        "-gf",
        "--glmPCAfamily",
        help="The choice of exponential family distribution to use for glmPCA method. (Default: %(default)s)",
    ),
    numberOfProcessors: str = typer.Option(
        "max",
        "-p",
        "--numberOfProcessors",
        callback=normalize_processors,
        help='Number of processors to use. Type "max/2" to use half the maximum number of processors or "max" to use all available processors. (Default: "max")',
    ),
    verbose: bool = typer.Option(False, "-v", "--verbose", help="Set to see processing messages."),
) -> int:
    _common(input, outFile, nComps, nNeighbors, numberOfProcessors, verbose)
    log_parameters(glmPCAfamily=glmPCAfamily)
    return 0


if __name__ == "__main__":
    app()
