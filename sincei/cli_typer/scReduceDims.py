from __future__ import annotations

import typer

from ._common_args import INPUT_OUTPUT_OPTS, OTHER_OPTS, GLMPCAFamily, log_parameters, override, preprocess_args

DESCRIPTION = (
    "Perform dimensionality reduction and UMAP projection on a cell-by-feature matrix.\n\n"
    "``scReduceDims`` performs dimensionality reduction on the input count matrix (output of scCountReads)"
    "and 2D projection (UMAP) of the cells. The result is an updated h5ad object with the dimensionality"
    "reduction results and UMAP coordinates in the ``.obsm`` field.\n\n"
    "``scReduceDims`` provides the following dimensionality reduction methods:\n"
    "* LSA: Latent Semantic Analysis.\n"
    "* LDA: Latent Dirichlet Allocation.\n"
    "* logPCA: Principal Component Analysis preceded by a logarithm transform.\n"
    "* glmPCA: generalized PCA, with an exponential family distribution such as Poisson, Bernoulli, etc."
)


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)

lsa_app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help="Use Latent Semantic Analysis (LSA).",
    context_settings={"help_option_names": []},
)
lda_app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help="Use Latent Dirichlet Allocation (LDA).",
    context_settings={"help_option_names": []},
)
logpca_app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help="Use log transform + PCA.",
    context_settings={"help_option_names": []},
)
glmpca_app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help="Use generalized PCA.",
    context_settings={"help_option_names": []},
)

app.add_typer(lsa_app, name="LSA")
app.add_typer(lda_app, name="LDA")
app.add_typer(logpca_app, name="logPCA")
app.add_typer(glmpca_app, name="glmPCA")

_REDUCTION = "Dimensionality reduction options"

# Options shared by the dimensionality-reduction subcommands.
N_COMPS = typer.Option(
    20,
    "-n",
    "--nComps",
    rich_help_panel=_REDUCTION,
    help="Number of principal components or topics to reduce the dimensionality to. Use a higher number for samples "
    "with more expected heterogeneity.",
)
N_NEIGHBORS = typer.Option(
    30,
    "-nk",
    "--nNeighbors",
    rich_help_panel=_REDUCTION,
    help="Number of nearest neighbours to consider for UMAP. Choose this considering the total number of cells and "
    "the expected number of clusters; smaller numbers lead to more fragmented clusters.",
)
BINARIZE = typer.Option(
    False,
    "--binarize",
    rich_help_panel=_REDUCTION,
    help="Binarize the counts per region before dimensionality reduction.",
)
OUT_FILE_TRAINED_MODEL = typer.Option(
    None,
    "-om",
    "--outFileTrainedModel",
    rich_help_panel=_REDUCTION,
    help="The output file for the trained model. The saved model can be used later to embed/compare new cells to the "
    "existing cluster of cells.",
)


@app.callback(invoke_without_command=True)
def main(ctx: typer.Context, help: bool = override(OTHER_OPTS["help"], rich_help_panel="Options")) -> int:
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
        raise typer.Exit()
    return 0


def _common(**parameters: object) -> None:
    log_parameters(**parameters)


@lsa_app.callback(invoke_without_command=True)
def lsa(
    # Input / Output options
    input: str = INPUT_OUTPUT_OPTS["h5ad_file"],
    out_file: str = INPUT_OUTPUT_OPTS["out_file"],
    # Dimensionality reduction options
    n_comps: int = N_COMPS,
    n_neighbors: int = N_NEIGHBORS,
    binarize: bool = BINARIZE,
    out_file_trained_model: str = OUT_FILE_TRAINED_MODEL,
    # Other options
    number_of_processors: str = OTHER_OPTS["number_of_processors"],
    verbose: bool = OTHER_OPTS["verbose"],
    help: bool = OTHER_OPTS["help"],
) -> int:
    _common(
        input=input,
        out_file=out_file,
        n_comps=n_comps,
        n_neighbors=n_neighbors,
        number_of_processors=number_of_processors,
        verbose=verbose,
    )
    log_parameters(binarize=binarize, out_file_trained_model=out_file_trained_model)
    return 0


@lda_app.callback(invoke_without_command=True)
def lda(
    # Input / Output options
    input: str = INPUT_OUTPUT_OPTS["h5ad_file"],
    out_file: str = INPUT_OUTPUT_OPTS["out_file"],
    # Dimensionality reduction options
    n_comps: int = N_COMPS,
    n_neighbors: int = N_NEIGHBORS,
    binarize: bool = BINARIZE,
    out_file_trained_model: str = OUT_FILE_TRAINED_MODEL,
    # LDA options
    n_passes: int = typer.Option(
        5,
        "--nPasses",
        rich_help_panel="LDA options",
        help="Number of passes through the corpus for LDA model fitting.",
    ),
    n_iterations: int = typer.Option(
        50,
        "--nIterations",
        rich_help_panel="LDA options",
        help="Number of iterations per pass for LDA model fitting.",
    ),
    alpha: float = typer.Option(
        1.0,
        "--alpha",
        rich_help_panel="LDA options",
        help="Prior to initialize cell-topic vectors.",
    ),
    eta: float = typer.Option(
        0.1,
        "--eta",
        rich_help_panel="LDA options",
        help="Prior to initialize feature-topic vectors.",
    ),
    gamma_threshold: float = typer.Option(
        None,
        "--gammaThreshold",
        rich_help_panel="LDA options",
        help="Minimum change in the topic matrix to stop the LDA model fitting. If not given, the model is fit for "
        "the number of passes and iterations specified above.",
    ),
    # Other options
    number_of_processors: str = OTHER_OPTS["number_of_processors"],
    verbose: bool = OTHER_OPTS["verbose"],
    help: bool = OTHER_OPTS["help"],
) -> int:
    _common(
        input=input,
        out_file=out_file,
        n_comps=n_comps,
        n_neighbors=n_neighbors,
        number_of_processors=number_of_processors,
        verbose=verbose,
    )
    log_parameters(
        binarize=binarize,
        out_file_trained_model=out_file_trained_model,
        n_passes=n_passes,
        n_iterations=n_iterations,
        alpha=alpha,
        eta=eta,
        gamma_threshold=gamma_threshold,
    )
    return 0


@logpca_app.callback(invoke_without_command=True)
def logpca(
    # Input / Output options
    input: str = INPUT_OUTPUT_OPTS["h5ad_file"],
    out_file: str = INPUT_OUTPUT_OPTS["out_file"],
    # Dimensionality reduction options
    n_comps: int = N_COMPS,
    n_neighbors: int = N_NEIGHBORS,
    # Other options
    number_of_processors: str = OTHER_OPTS["number_of_processors"],
    verbose: bool = OTHER_OPTS["verbose"],
    help: bool = OTHER_OPTS["help"],
) -> int:
    _common(
        input=input,
        out_file=out_file,
        n_comps=n_comps,
        n_neighbors=n_neighbors,
        number_of_processors=number_of_processors,
        verbose=verbose,
    )
    return 0


@glmpca_app.callback(invoke_without_command=True)
def glmpca(
    # Input / Output options
    input: str = INPUT_OUTPUT_OPTS["h5ad_file"],
    out_file: str = INPUT_OUTPUT_OPTS["out_file"],
    # Dimensionality reduction options
    n_comps: int = N_COMPS,
    n_neighbors: int = N_NEIGHBORS,
    # glmPCA options
    glmpca_family: GLMPCAFamily = typer.Option(
        GLMPCAFamily.poisson,
        "-gf",
        "--glmPCAfamily",
        metavar="FAMILY",
        rich_help_panel="glmPCA options",
        help="The choice of exponential family distribution to use for the glmPCA method.\n\n"
        "One of: [bold yellow]poisson[/bold yellow], [bold yellow]nb[/bold yellow], "
        "[bold yellow]mult[/bold yellow], [bold yellow]bern[/bold yellow].",
    ),
    # Other options
    number_of_processors: str = OTHER_OPTS["number_of_processors"],
    verbose: bool = OTHER_OPTS["verbose"],
    help: bool = OTHER_OPTS["help"],
) -> int:
    _common(
        input=input,
        out_file=out_file,
        n_comps=n_comps,
        n_neighbors=n_neighbors,
        number_of_processors=number_of_processors,
        verbose=verbose,
    )
    log_parameters(glmpca_family=glmpca_family)
    return 0


def cli() -> None:
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
