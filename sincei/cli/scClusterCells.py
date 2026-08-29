from __future__ import annotations

from typing import Annotated

import typer

from ._common_args import (
    AVAILABLE_PROCESSORS,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    PLOT_OPTS,
    DimRed,
    GLMPCAFamily,
    PlotFileFormat,
    configure_logging,
    log_parameters,
    preprocess_args,
)

DESCRIPTION = (
    "Cluster cells from a cell-by-feature matrix.\n\n"
    "``scClusterCells`` clusters cells based on the input count matrix (output of "
    "``scCountReads``) and performs dimensionality reduction, community detection and "
    "2D projection (UMAP) of the cells. The result is an updated h5ad object, and "
    "(optionally) a plot file and a .tsv file with UMAP coordinates and the "
    "corresponding cluster id for each barcode.\n\n"
    "If the input AnnData object already contains a dimensionality reduction of the "
    "same type as requested, the existing reduction is used instead of computing a "
    "new one.\n\n"
    "``scClusterCells`` provides the following dimensionality reduction methods:\n"
    "* LSA: Latent Semantic Analysis.\n"
    "* LDA: Latent Dirichlet Allocation.\n"
    "* logPCA: Principal Component Analysis preceded by a logarithm transform.\n"
    "* glmPCA: generalized PCA, with an exponential family distribution such as "
    "Poisson, Bernoulli, etc."
)


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)

_REDUCTION = "Dimensionality reduction options"
_LDA = "LDA options"
_GLMPCA = "glmPCA options"
_CLUSTERING = "Clustering options"


@app.callback(invoke_without_command=True)
def main(
    # Input / Output options
    input: Annotated[str, INPUT_OUTPUT_OPTS["h5ad_file"]],
    out_file: Annotated[str, INPUT_OUTPUT_OPTS["out_file"]],
    # Dimensionality reduction options
    method: Annotated[
        DimRed,
        typer.Option(
            "-m",
            "--method",
            metavar="METHOD",
            rich_help_panel=_REDUCTION,
            help=(
                "The dimensionality reduction method to use before clustering cells."
                "\n\n"
                "One of: "
                "[bold yellow]LSA[/bold yellow], "
                "[bold yellow]LDA[/bold yellow], "
                "[bold yellow]logPCA[/bold yellow], "
                "[bold yellow]glmPCA[/bold yellow]."
            ),
        ),
    ] = DimRed.LSA,
    n_prin_comps: Annotated[
        int,
        typer.Option(
            "-n",
            "--nPrinComps",
            rich_help_panel=_REDUCTION,
            help=(
                "Number of principal components or topics to reduce the dimensionality "
                "to. Use a higher number for samples with more expected heterogeneity."
            ),
        ),
    ] = 20,
    n_neighbors: Annotated[
        int,
        typer.Option(
            "-nk",
            "--nNeighbors",
            rich_help_panel=_REDUCTION,
            help=(
                "Number of nearest neighbours to consider for clustering and UMAP. "
                "Choose this considering the total number of cells and the expected "
                "number of clusters; smaller numbers lead to more fragmented clusters."
            ),
        ),
    ] = 30,
    binarize: Annotated[
        bool,
        typer.Option(
            "--binarize",
            rich_help_panel=_REDUCTION,
            help=(
                "Binarize the counts per region before dimensionality reduction. Only "
                "used with [bold yellow]LSA[/bold yellow] and "
                "[bold yellow]LDA[/bold yellow]."
            ),
        ),
    ] = False,
    out_file_trained_model: Annotated[
        str | None,
        typer.Option(
            "-om",
            "--outFileTrainedModel",
            rich_help_panel=_REDUCTION,
            help=(
                "The output file for the trained model. The saved model can be used "
                "later to embed/compare new cells to the existing cluster of cells."
            ),
        ),
    ] = None,
    recompute_reduction: Annotated[
        bool,
        typer.Option(
            "--recomputeReduction",
            rich_help_panel=_REDUCTION,
            help=(
                "Recompute the dimensionality reduction even if a precomputed version "
                "exists."
            ),
        ),
    ] = False,
    # LDA options
    n_passes: Annotated[
        int,
        typer.Option(
            "--nPasses",
            rich_help_panel=_LDA,
            help=(
                "Number of passes through the corpus for LDA model fitting. Only used "
                "with [bold yellow]LDA[/bold yellow]."
            ),
        ),
    ] = 5,
    n_iterations: Annotated[
        int,
        typer.Option(
            "--nIterations",
            rich_help_panel=_LDA,
            help=(
                "Number of iterations per pass for LDA model fitting. Only used with "
                "[bold yellow]LDA[/bold yellow]."
            ),
        ),
    ] = 50,
    alpha: Annotated[
        float,
        typer.Option(
            "--alpha",
            rich_help_panel=_LDA,
            help=(
                "Prior to initialize cell-topic vectors. Only used with "
                "[bold yellow]LDA[/bold yellow]."
            ),
        ),
    ] = 1.0,
    eta: Annotated[
        float,
        typer.Option(
            "--eta",
            rich_help_panel=_LDA,
            help=(
                "Prior to initialize feature-topic vectors. Only used with "
                "[bold yellow]LDA[/bold yellow]."
            ),
        ),
    ] = 0.1,
    gamma_threshold: Annotated[
        float | None,
        typer.Option(
            "--gammaThreshold",
            rich_help_panel=_LDA,
            help=(
                "Minimum change in the topic matrix to stop the LDA model fitting. If "
                "not given, the model is fit for the number of passes and iterations "
                "specified above. Only used with [bold yellow]LDA[/bold yellow]."
            ),
        ),
    ] = None,
    # glmPCA options
    glmpca_family: Annotated[
        GLMPCAFamily,
        typer.Option(
            "-gf",
            "--glmPCAfamily",
            metavar="FAMILY",
            rich_help_panel=_GLMPCA,
            help=(
                "The choice of exponential family distribution to use for the glmPCA "
                "method. Only used with [bold yellow]glmPCA[/bold yellow].\n\n"
                "One of: [bold yellow]poisson[/bold yellow], "
                "[bold yellow]nb[/bold yellow], [bold yellow]mult[/bold yellow], "
                "[bold yellow]bern[/bold yellow]."
            ),
        ),
    ] = GLMPCAFamily.poisson,
    # Clustering options
    out_file_umap: Annotated[
        str | None,
        typer.Option(
            "-op",
            "--outFileUMAP",
            rich_help_panel=_CLUSTERING,
            help=(
                "The output plot file (for UMAP). If specified, a 4-column .tsv file "
                "with the same prefix is also created with the cell IDs, raw UMAP "
                "coordinates (UMAP1 and UMAP2) and Leiden cluster number."
            ),
        ),
    ] = None,
    cluster_resolution: Annotated[
        float,
        typer.Option(
            "-cr",
            "--clusterResolution",
            rich_help_panel=_CLUSTERING,
            help=(
                "Resolution parameter for Leiden clustering. Values lower than 1.0 "
                "result in fewer clusters, while higher values lead to splitting of "
                "clusters. In most cases the optimum is between 0.8 and 1.2."
            ),
        ),
    ] = 1.0,
    # Plot options
    plot_width: Annotated[float, PLOT_OPTS["plot_width"]] = 10.0,
    plot_height: Annotated[float, PLOT_OPTS["plot_height"]] = 10.0,
    plot_file_format: Annotated[
        PlotFileFormat, PLOT_OPTS["plot_file_format"]
    ] = PlotFileFormat.png,
    # Other options
    number_of_processors: Annotated[
        int, OTHER_OPTS["number_of_processors"]
    ] = AVAILABLE_PROCESSORS,
    verbose: Annotated[bool, OTHER_OPTS["verbose"]] = False,
    help: Annotated[bool, OTHER_OPTS["help"]] = False,
) -> int:
    log_parameters(
        input=input,
        out_file=out_file,
        method=method,
        n_prin_comps=n_prin_comps,
        n_neighbors=n_neighbors,
        binarize=binarize,
        recompute_reduction=recompute_reduction,
        out_file_trained_model=out_file_trained_model,
        n_passes=n_passes,
        n_iterations=n_iterations,
        alpha=alpha,
        eta=eta,
        gamma_threshold=gamma_threshold,
        glmpca_family=glmpca_family,
        out_file_umap=out_file_umap,
        cluster_resolution=cluster_resolution,
        plot_width=plot_width,
        plot_height=plot_height,
        plot_file_format=plot_file_format,
        number_of_processors=number_of_processors,
        verbose=verbose,
    )
    return 0


def cli() -> None:
    configure_logging()
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
