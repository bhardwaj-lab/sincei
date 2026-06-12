from __future__ import annotations

import typer

from ._common_args import (
    BAM_OPTS,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    log_parameters,
    override,
    preprocess_args,
)

DESCRIPTION = (
    "Merge AnnData files from different modalities into a MuData object.\n\n"
    "``scCombineMods`` combines multiple count matrices (output of scCountReads) of different data "
    "modalities (e.g. gene expression, chromatin accessibility, histone modifications) into one. "
    "The result is a .h5mu (MuData) file containing each of the data modalities provided.\n"
    "*NOTE*: this does not perform any 'batch effect correction' or 'integration' of data from different "
    "technologies."
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
    input: list[str] = INPUT_OUTPUT_OPTS["h5ad_files"],
    labels: list[str] = override(BAM_OPTS["labels"], default=...),
    out_file: str = INPUT_OUTPUT_OPTS["out_file"],
    # Other options
    number_of_processors: str = OTHER_OPTS["number_of_processors"],
    verbose: bool = OTHER_OPTS["verbose"],
    help: bool = OTHER_OPTS["help"],
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
