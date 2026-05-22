from __future__ import annotations

import typer

from .shared import log_parameters, version_option

DESCRIPTION = """
Merge AnnData files from different modalities into a MuData object.
"""


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help=DESCRIPTION,
    context_settings={"help_option_names": ["-h", "--help"]},
)


@version_option("scCombineMods")
@app.callback(invoke_without_command=True)
def main(
    input: list[str] = typer.Option(
        ..., "-i", "--input", metavar=".h5ad", help="List of sincei-generated input .h5ad files separated by spaces."
    ),
    labels: list[str] = typer.Option(
        ...,
        "-l",
        "--labels",
        metavar="sample",
        help="User defined labels instead of default labels from file names. Multiple labels have to be separated by a space, e.g. ``--labels sample1 sample2 sample3``.",
    ),
    outFile: str = typer.Option(
        ...,
        "-o",
        "--outFile",
        help="The file to write results to. For `scFilterStats`, `scFilterBarcodes` and `scJSD`, the output file is a .tsv file. For other tools, the output file is an updated .h5ad or .h5mu file with the result of the requested operation.",
    ),
    numberOfProcessors: str = typer.Option(
        "max",
        "-p",
        "--numberOfProcessors",
        help='Number of processors to use. Type "max/2" to use half the maximum number of processors or "max" to use all available processors. (Default: "max")',
    ),
    verbose: bool = typer.Option(False, "-v", "--verbose", help="Set to see processing messages."),
) -> int:
    log_parameters(input=input, labels=labels, outFile=outFile, numberOfProcessors=numberOfProcessors, verbose=verbose)
    return 0


if __name__ == "__main__":
    app()
