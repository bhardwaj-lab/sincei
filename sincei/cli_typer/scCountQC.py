from __future__ import annotations

from typing import Annotated

import typer

from ._common_args import (
    AVAILABLE_PROCESSORS,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    log_parameters,
    preprocess_args,
)

DESCRIPTION = (
    "Perform quality control and filter cells and regions from a cell-by-feature "
    "matrix.\n\n"
    "``scCountQC`` calculates multiple quality controls metrics on the input .h5ad "
    "file (output of scCountReads) and (optionally) filters the input file based on "
    "filterCellArgs/filterRegionArgs. The output is either an updated .h5ad object (if "
    "filtering is requested) or the filtering metrics (--outMetrics)."
)


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)

_QC = "QC options"


@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    # Input and output options
    input: Annotated[str, INPUT_OUTPUT_OPTS["h5ad_file"]],
    out_file: Annotated[str, INPUT_OUTPUT_OPTS["out_file"]],
    # QC options
    describe: Annotated[
        bool,
        typer.Option(
            "-d",
            "--describe",
            rich_help_panel=_QC,
            help="Print a list of cell and region metrics available for QC/filtering.",
        ),
    ] = False,
    out_metrics: Annotated[
        str | None,
        typer.Option(
            "-om",
            "--outMetrics",
            rich_help_panel=_QC,
            help=(
                "Prefix of the output file with calculated QC metrics. If given, the "
                "cell metrics are printed in <prefix>.cell.tsv and region metrics as "
                "<prefix>.region.tsv."
            ),
        ),
    ] = None,
    filter_cell_args: Annotated[
        str | None,
        typer.Option(
            "-fc",
            "--filterCellArgs",
            rich_help_panel=_QC,
            help=(
                'List of arguments to filter cells. The format is "arg_name: minvalue, '
                'maxvalue; arg_name: minvalue, maxvalue; ...." where arg_name is a '
                'cell QC metric present in the input h5ad file. Run with "--describe" '
                "to view all available metrics. The two values are used as lower and "
                "upper bounds to filter cells."
            ),
        ),
    ] = None,
    filter_region_args: Annotated[
        str | None,
        typer.Option(
            "-fr",
            "--filterRegionArgs",
            rich_help_panel=_QC,
            help=(
                'List of arguments to filter regions. The format is "arg_name: '
                'minvalue, maxvalue; arg_name: minvalue; ...." where arg_name is a '
                "region QC metric present in the input h5ad file. Run with "
                '"--describe" to view all available metrics. The two values are used '
                "as lower and upper bounds to filter regions."
            ),
        ),
    ] = None,
    region_blacklist: Annotated[
        list[str] | None,
        typer.Option(
            "-rb",
            "--regionBlacklist",
            metavar="BED",
            rich_help_panel=_QC,
            help=(
                "A BED or GTF file containing regions that should be excluded from all "
                "analyses. Regions in the anndata object that overlap with blacklisted "
                "regions will be removed."
            ),
        ),
    ] = None,
    cell_blacklist: Annotated[
        str | None,
        typer.Option(
            "-cb",
            "--cellBlacklist",
            rich_help_panel=_QC,
            help=(
                "A list of barcodes to be excluded from the clustering. The barcodes "
                "(along with sample labels) must be present in the input object."
            ),
        ),
    ] = None,
    chrom_blacklist: Annotated[
        list[str] | None,
        typer.Option(
            "-chb",
            "--chromBlacklist",
            metavar="CHR",
            rich_help_panel=_QC,
            help="A space separated list of chromosomes to exclude, e.g. chrM chrUn.",
        ),
    ] = None,
    # Other options
    number_of_processors: Annotated[
        int, OTHER_OPTS["number_of_processors"]
    ] = AVAILABLE_PROCESSORS,
    verbose: Annotated[bool, OTHER_OPTS["verbose"]] = False,
    help: Annotated[bool, OTHER_OPTS["help"]] = False,
) -> int:
    if ctx.invoked_subcommand is None:
        log_parameters(
            input=input,
            out_file=out_file,
            describe=describe,
            out_metrics=out_metrics,
            filter_cell_args=filter_cell_args,
            filter_region_args=filter_region_args,
            region_blacklist=region_blacklist,
            cell_blacklist=cell_blacklist,
            chrom_blacklist=chrom_blacklist,
            number_of_processors=number_of_processors,
            verbose=verbose,
        )
    return 0


def cli() -> None:
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
