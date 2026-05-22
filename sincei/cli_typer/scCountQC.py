from __future__ import annotations

import typer

from .shared import log_parameters, version_option

DESCRIPTION = """
Perform quality control and filter cells and regions from a cell-by-feature matrix.
"""


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help=DESCRIPTION,
    context_settings={"help_option_names": ["-h", "--help"]},
)


@version_option("scCountQC")
@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    input: str = typer.Option(
        ..., "-i", "--input", metavar=".h5ad", help="sincei-generated input file in .h5ad format."
    ),
    outFile: str = typer.Option(
        ...,
        "-o",
        "--outFile",
        help="The file to write results to. For `scFilterStats`, `scFilterBarcodes` and `scJSD`, the output file is a .tsv file. For other tools, the output file is an updated .h5ad or .h5mu file with the result of the requested operation.",
    ),
    describe: bool = typer.Option(
        False, "-d", "--describe", help="Print a list of cell and region metrics available for QC/filtering."
    ),
    outMetrics: str = typer.Option(
        None,
        "-om",
        "--outMetrics",
        help="Prefix of the output file with calculated QC metrics. If given, the cell metrics are printed in <prefix>.cell.tsv and region metrics as <prefix>.region.tsv",
    ),
    filterCellArgs: str = typer.Option(
        None,
        "-fc",
        "--filterCellArgs",
        help='List of arguments to filter cells. The format is "arg_name: minvalue, maxvalue; arg_name: minvalue, maxvalue; ...." where arg_name is the QC metric for cells present in the input h5ad file. In order to view all available cell filtering metrics, run scCountFilter with "--describe". The two arguments are supplied (minvalue, maxvalue) they are used as lower and upper bounds to filter cells. Make sure they are float/integer numbers.',
    ),
    filterRegionArgs: str = typer.Option(
        None,
        "-fr",
        "--filterRegionArgs",
        help='List of arguments to filter regions. The format is "arg_name: minvalue, maxvalue; arg_name: minvalue; ...." where arg_name is the QC metric for regions present in the input h5ad file. In order to view all available cell filtering metrics, run scCountFilter with "--describe". The two arguments are supplied (minvalue, maxvalue) they are used as lower and upper bounds to filter cells. Make sure they are float/integer numbers.',
    ),
    region_blacklist: list[str] = typer.Option(
        None,
        "-rb",
        "--region_blacklist",
        metavar="BED",
        help="A BED or GTF file containing regions that should be excluded from all analyses. Regions in the anndata object that overlap with blacklisted regions will be removed.",
    ),
    cell_blacklist: str = typer.Option(
        None,
        "-cb",
        "--cell_blacklist",
        help="A list of barcodes to be excluded for the clustering. The barcodes (along with sample labels) must be present in the input object.",
    ),
    chrom_blacklist: list[str] = typer.Option(
        None,
        "-chb",
        "--chrom_blacklist",
        metavar="CHR1",
        help="A space separated list of chromosomes to exclude. eg. chrM chrUn",
    ),
    numberOfProcessors: str = typer.Option(
        "max",
        "-p",
        "--numberOfProcessors",
        help='Number of processors to use. Type "max/2" to use half the maximum number of processors or "max" to use all available processors. (Default: "max")',
    ),
    verbose: bool = typer.Option(False, "-v", "--verbose", help="Set to see processing messages."),
) -> int:
    if ctx.invoked_subcommand is None:
        log_parameters(
            input=input,
            outFile=outFile,
            describe=describe,
            outMetrics=outMetrics,
            filterCellArgs=filterCellArgs,
            filterRegionArgs=filterRegionArgs,
            region_blacklist=region_blacklist,
            cell_blacklist=cell_blacklist,
            chrom_blacklist=chrom_blacklist,
            numberOfProcessors=numberOfProcessors,
            verbose=verbose,
        )
    return 0


if __name__ == "__main__":
    app()
