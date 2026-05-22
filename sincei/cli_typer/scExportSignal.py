from __future__ import annotations

from typing import Annotated

import typer

from .shared import log_parameters, normalize_processors, normalize_region, version_option

DESCRIPTION = """
scExportSignal exports signal from a .h5ad object in a few output formats.
"""


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help=DESCRIPTION,
    context_settings={"help_option_names": ["-h", "--help"]},
)


@version_option("scExportSignal")
@app.callback(invoke_without_command=True)
def main(
    input: Annotated[
        str, typer.Option("-i", "--input", metavar=".h5ad", help="sincei-generated input file in .h5ad format.")
    ],
    outFile: Annotated[str, typer.Option("-o", "--outFile", help="The file to write results to.")],
    outFileFormat: Annotated[
        str,
        typer.Option(
            "--outFileFormat",
            help="Output file format for `scExportSignal`. "
            '"bm" refers to the BedGraphMatrix format, useful for single-cell data visualization with pyGenomeTracks. It stores data in dense format, so we recommend choosing a region to export instead of the whole dataset. '
            '"mtx" refers to the MatrixMarket sparse-matrix format. The output in this case would be <prefix>.counts.mtx, along with <prefix>.rownames.txt and <prefix>.colnames.txt'
            '"loom" refers to the loom file format, an hddf5-based legacy format for single-cell genomics data.',
        ),
    ],
    region: str = typer.Option(None, "-r", "--region", callback=normalize_region),
    numberOfProcessors: str = typer.Option(
        "max",
        "-p",
        "--numberOfProcessors",
        callback=normalize_processors,
        help='Number of processors to use. Type "max/2" to use half the maximum number of processors or "max" to use all available processors. (Default: "max")',
    ),
    verbose: bool = typer.Option(False, "-v", "--verbose", help="Set to see processing messages."),
) -> int:
    log_parameters(
        input=input,
        outFile=outFile,
        region=region,
        outFileFormat=outFileFormat,
        numberOfProcessors=numberOfProcessors,
        verbose=verbose,
    )
    return 0


if __name__ == "__main__":
    app()
