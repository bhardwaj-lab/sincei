from __future__ import annotations

from typing import Annotated

import typer

from ._common_args import (
    AVAILABLE_PROCESSORS,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    PLOT_OPTS,
    DimRed,
    PlotFileFormat,
    log_parameters,
    preprocess_args,
)

DESCRIPTION = (
    "Cluster cells from a cell-by-feature matrix using the Leiden algorithm.\n\n"
    "``scClusterCells`` clusters cells based on the dimensionality reduction in the "
    "input h5ad file. The result is an updated h5ad object, and (optionally) a plot "
    "and a .tsv file with UMAP coordinates and corresponding cluster id for each cell."
)


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)

_CLUSTERING = "Clustering options"


@app.callback(invoke_without_command=True)
def main(
    # Input / Output options
    input: Annotated[str, INPUT_OUTPUT_OPTS["h5ad_file"]],
    out_file: Annotated[str, INPUT_OUTPUT_OPTS["out_file"]],
    # Clustering options
    out_file_umap: Annotated[
        str | None,
        typer.Option(
            "-ou",
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
    dim_red: Annotated[
        DimRed | None,
        typer.Option(
            "--dimRed",
            metavar="METHOD",
            rich_help_panel=_CLUSTERING,
            help=(
                "Dimensionality reduction modality to perform Leiden clustering on. If "
                "not given, the program searches the ``.obsm`` field of the input h5ad "
                "for the output of ``scReduceDims`` in order of preference: LSA, LDA, "
                "logPCA, glmPCA.\n\n"
                "One of: [bold yellow]LSA[/bold yellow], [bold yellow]LDA[/bold "
                "yellow], [bold yellow]logPCA[/bold yellow], "
                "[bold yellow]glmPCA[/bold yellow]."
            ),
        ),
    ] = None,
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
        out_file_umap=out_file_umap,
        cluster_resolution=cluster_resolution,
        dim_red=dim_red,
        plot_width=plot_width,
        plot_height=plot_height,
        plot_file_format=plot_file_format,
        number_of_processors=number_of_processors,
        verbose=verbose,
    )
    return 0


def cli() -> None:
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
