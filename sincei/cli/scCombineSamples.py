from __future__ import annotations

from typing import Annotated

import typer

from ._common_args import (
    _IO,
    AVAILABLE_PROCESSORS,
    BAM_OPTS,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    log_parameters,
    override,
    preprocess_args,
)

DESCRIPTION = (
    "Concatenate/merge AnnData files from different samples.\n\n"
    "``scCombineSamples`` combines multiple count matrices (output of scCountReads) "
    "into one. Only features present in all matrices will be kept. The result is a "
    ".h5ad (AnnData) file containing the combined count matrix.\n"
    "*NOTE*: it doesn't perform any 'batch effect correction' or 'integration' of "
    "data from different technologies."
)


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
    input: Annotated[list[str], INPUT_OUTPUT_OPTS["h5ad_files"]],
    out_file: Annotated[str, INPUT_OUTPUT_OPTS["out_file"]],
    labels: Annotated[
        list[str] | None, override(BAM_OPTS["labels"], rich_help_panel=_IO)
    ] = None,
    smart_labels: Annotated[
        bool, override(BAM_OPTS["smart_labels"], rich_help_panel=_IO)
    ] = False,
    # Other options
    number_of_processors: Annotated[
        int, OTHER_OPTS["number_of_processors"]
    ] = AVAILABLE_PROCESSORS,
    verbose: Annotated[bool, OTHER_OPTS["verbose"]] = False,
    help: Annotated[bool, OTHER_OPTS["help"]] = False,
) -> int:
    log_parameters(
        input=input,
        labels=labels,
        out_file=out_file,
        smart_labels=smart_labels,
        number_of_processors=number_of_processors,
        verbose=verbose,
    )
    return 0


def cli() -> None:
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
