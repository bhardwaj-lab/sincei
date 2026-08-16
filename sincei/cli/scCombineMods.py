from __future__ import annotations

from typing import Annotated

import typer

from ._common_args import (
    AVAILABLE_PROCESSORS,
    BAM_OPTS,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    log_parameters,
    preprocess_args,
)

DESCRIPTION = (
    "Merge AnnData files from different modalities into a MuData object.\n\n"
    "``scCombineMods`` combines multiple count matrices (output of scCountReads) of "
    "different data modalities (e.g. gene expression, chromatin accessibility, "
    "histone modifications) into one. The result is a .h5mu (MuData) file containing "
    "each of the data modalities provided.\n"
    "*NOTE*: this does not perform any 'batch effect correction' or 'integration' of "
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
    labels: Annotated[list[str], BAM_OPTS["labels"]],
    out_file: Annotated[str, INPUT_OUTPUT_OPTS["out_file"]],
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
        number_of_processors=number_of_processors,
        verbose=verbose,
    )
    return 0


def cli() -> None:
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
