from __future__ import annotations

from typing import Annotated

import typer

from ._common_args import (
    _IO,
    AVAILABLE_PROCESSORS,
    BAM_OPTS,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    CombineMethod,
    configure_logging,
    log_parameters,
    override,
    preprocess_args,
)

DESCRIPTION = (
    "Combine multiple count matrices into one.\n\n"
    "``scCombineCounts`` combines multiple count matrices (output of "
    "``scCountReads``) into one, either assuming they are different samples "
    "([bold yellow]multi-sample[/bold yellow]) or different measurements on the same "
    "set of cells ([bold yellow]multi-modal[/bold yellow]).\n\n"
    "* [bold yellow]multi-sample[/bold yellow]: each sample is independent but was "
    "counted on the same features, so the tool looks for feature overlaps and not "
    "for barcode overlaps. Only features present in every matrix are kept. The "
    "result is a .h5ad (AnnData) file.\n"
    "* [bold yellow]multi-modal[/bold yellow]: the counts come from different assays "
    "on the same set of cells, so the tool looks for cell-barcode overlaps and not "
    "for feature overlaps. The result is a .h5mu (MuData) file, which "
    "``scClusterCells`` can use for multi-modal clustering.\n\n"
    "*NOTE*: this performs no 'batch effect correction' or 'integration' of data "
    "from different technologies, which needs more sophisticated methods."
)


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)

_GENERAL = "General options"


@app.callback(invoke_without_command=True)
def main(
    # Input / Output options
    input: Annotated[list[str], INPUT_OUTPUT_OPTS["h5ad_files"]],
    out_file: Annotated[str, INPUT_OUTPUT_OPTS["out_file"]],
    labels: Annotated[
        list[str] | None, override(BAM_OPTS["labels"], rich_help_panel=_IO)
    ] = None,
    smart_labels: Annotated[
        bool,
        override(
            BAM_OPTS["smart_labels"],
            rich_help_panel=_IO,
            # The shared text says "BAM files"; this tool takes .h5ad files.
            help=(
                "Instead of manually specifying labels for the input files, use the "
                "file name after removing the path and extension."
            ),
        ),
    ] = False,
    # General options
    method: Annotated[
        CombineMethod,
        typer.Option(
            "-m",
            "--method",
            metavar="METHOD",
            rich_help_panel=_GENERAL,
            help=(
                "How to merge the counts from the provided matrices.\n\n"
                "One of: [bold yellow]multi-sample[/bold yellow], "
                "[bold yellow]multi-modal[/bold yellow]."
            ),
        ),
    ] = CombineMethod.multi_sample,
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
        labels=labels,
        smart_labels=smart_labels,
        method=method,
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
