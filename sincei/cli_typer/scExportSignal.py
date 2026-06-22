from __future__ import annotations

import typer

from ._common_args import (
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    ExportFormat,
    log_parameters,
    preprocess_args,
)

DESCRIPTION = "Export .h5ad objects to other formats."


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)


@app.callback(invoke_without_command=True)
def main(
    # Input / Output options
    input: str = INPUT_OUTPUT_OPTS["h5ad_file"],
    out_file: str = INPUT_OUTPUT_OPTS["out_prefix"],
    out_file_format: ExportFormat = typer.Option(
        ...,
        "--outFileFormat",
        metavar="FORMAT",
        rich_help_panel="Export options",
        help="Output file format.\n\n"
        "[bold yellow]bm[/bold yellow]: BedGraphMatrix format, useful for single-cell visualization with pyGenomeTracks; "
        "stores data densely, so prefer exporting a region rather than the whole dataset.\n\n"
        "[bold yellow]mtx[/bold yellow]: MatrixMarket sparse format (<prefix>.counts.mtx plus <prefix>.rownames.txt and "
        "<prefix>.colnames.txt).\n\n"
        "[bold yellow]loom[/bold yellow]: loom hdf5-based legacy format for single-cell genomics data.",
    ),
    region: str | None = INPUT_OUTPUT_OPTS["region"],
    # Other options
    number_of_processors: str = OTHER_OPTS["number_of_processors"],
    verbose: bool = OTHER_OPTS["verbose"],
    help: bool = OTHER_OPTS["help"],
) -> int:
    log_parameters(
        input=input,
        out_file=out_file,
        out_file_format=out_file_format,
        region=region,
        number_of_processors=number_of_processors,
        verbose=verbose,
    )
    return 0


def cli() -> None:
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
