from __future__ import annotations

import typer

from .shared import log_parameters, normalize_processors, normalize_region, version_option

DESCRIPTION = """
Call variable chromatin regions (VCRs) from binned chromatin data.
"""


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help=DESCRIPTION,
    context_settings={"help_option_names": ["-h", "--help"]},
)


@version_option("scFindVCRs")
@app.callback(invoke_without_command=True)
def main(
    input: str = typer.Option(
        ..., "-i", "--input", metavar=".h5ad", help="sincei-generated input file in .h5ad format."
    ),
    region: str = typer.Option(None, "-r", "--region", callback=normalize_region),
    binSize: int = typer.Option(
        ..., "-bs", "--binSize", metavar="INT", help="The size of the bins in the input Anndata object."
    ),
    maxRegionSize: int = typer.Option(
        None,
        "-mr",
        "--maxRegionSize",
        metavar="INT",
        help="The maximum region size to be considered, in base pairs. Larger regions may increase compute time. Defaults to 100 times the bin size.",
    ),
    nKernels: int = typer.Option(
        20,
        "-nk",
        "--nKernels",
        metavar="INT",
        help="The number of kernels to use for the score map. More kernels generally lead to a better segmentation, but increase the computational cost.",
    ),
    penalties: list[float] = typer.Option(
        [0.05, 0.1, 0.5],
        "-pen",
        "--penalties",
        metavar="INT",
        help='Penalty value for change-point detection. Higher values result in fewer segments. Multiple values can be provided (separated by space). Each penalty value will produce a separate set of regions within which can be seperated from the output BED file by filtering on the "score" column.',
    ),
    outFile: str = typer.Option(
        ...,
        "-o",
        "--outFile",
        metavar=".bed",
        help='Name of the output file (BED format) with genome segmentation result. The penalty threshold that defines each segment is saved in the "score" column of the BED file, and the BED file can be filtered based on this column to obtain non-overlapping segments.',
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
    log_parameters(
        input=input,
        region=region,
        binSize=binSize,
        maxRegionSize=maxRegionSize,
        nKernels=nKernels,
        penalties=penalties,
        outFile=outFile,
        numberOfProcessors=numberOfProcessors,
        verbose=verbose,
    )
    return 0


if __name__ == "__main__":
    app()
