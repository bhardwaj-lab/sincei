from __future__ import annotations

import typer

from .shared import log_parameters, version_option

DESCRIPTION = """
Cluster cells from a cell-by-feature matrix using the Leiden algorithm.
"""


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help=DESCRIPTION,
    context_settings={"help_option_names": ["-h", "--help"]},
)


@version_option("scClusterCells")
@app.callback(invoke_without_command=True)
def main(
    input: str = typer.Option(
        ..., "-i", "--input", metavar=".h5ad", help="sincei-generated input file in .h5ad format."
    ),
    outFile: str = typer.Option(..., "-o", "--outFile", help="The file to write results to."),
    outFileUMAP: str = typer.Option(
        None,
        "-op",
        "--outFileUMAP",
        help="The output plot file (for UMAP). If you specify this option, a 4-column .tsv file with the same prefix is also created with the cell IDs, raw UMAP coordinates (UMAP1 and UMAP2) and Leiden cluster number.",
    ),
    clusterResolution: float = typer.Option(
        1.0,
        "-cr",
        "--clusterResolution",
        help="Resolution parameter for Leiden clustering. Values lower than 1.0 would result in less clusters, while higher values lead to splitting of clusters. In most cases, the optimum value would be between 0.8 and 1.2. (Default: %(default)s)",
    ),
    dimRed: str = typer.Option(
        None,
        "--dimRed",
        help="Dimensionality reduction modality to perform Leiden clustering on. If not given, the program searches in the ``.obsm`` field of the input h5ad for the output of ``scReduceDims`` in order of preference: LSA, LDA, logPCA, glmPCA.",
    ),
    plotWidth: float = typer.Option(10.0, "--plotWidth", help="Output plot width (in cm). (Default: %(default)s)"),
    plotHeight: float = typer.Option(10.0, "--plotHeight", help="Output plot height (in cm). (Default: %(default)s)"),
    plotFileFormat: str = typer.Option(
        "png",
        "--plotFileFormat",
        help="Image format type. If given, this option will be used to save the plot file. (Default: %(default)s)",
    ),
    numberOfProcessors: str = typer.Option(
        "max",
        "-p",
        "--numberOfProcessors",
        help='Number of processors to use. Type "max/2" to use half the maximum number of processors or "max" to use all available processors. (Default: "max")',
    ),
    verbose: bool = typer.Option(False, "-v", "--verbose", help="Set to see processing messages."),
) -> int:
    log_parameters(
        input=input,
        outFile=outFile,
        outFileUMAP=outFileUMAP,
        clusterResolution=clusterResolution,
        dimRed=dimRed,
        plotWidth=plotWidth,
        plotHeight=plotHeight,
        plotFileFormat=plotFileFormat,
        numberOfProcessors=numberOfProcessors,
        verbose=verbose,
    )
    return 0


if __name__ == "__main__":
    app()
