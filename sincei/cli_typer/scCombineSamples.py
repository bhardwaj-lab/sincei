from __future__ import annotations

import typer

from ._common_args import (
    BAM_OPTS,
    INPUT_OUTPUT_OPTS,
    _IO,
    OTHER_OPTS,
    log_parameters,
    override,
    preprocess_args,
)

DESCRIPTION = (
    "Concatenate/merge AnnData files from different samples.\n\n"
    "``scCombineSamples`` combines multiple count matrices (output of scCountReads) into one. Only"
    "features present in all matrices will be kept. The result is a .h5ad (AnnData) file containing the"
    "combined count matrix.\n"
    "*NOTE*: it doesn't perform any 'batch effect correction' or 'integration' of data from different"
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
    out_file: str = INPUT_OUTPUT_OPTS["out_file"],
    labels: list[str] = override(BAM_OPTS["labels"], rich_help_panel=_IO),
    smart_labels: bool = override(BAM_OPTS["smart_labels"], rich_help_panel=_IO),
    # Other options
    number_of_processors: str = OTHER_OPTS["number_of_processors"],
    verbose: bool = OTHER_OPTS["verbose"],
    help: bool = OTHER_OPTS["help"],
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
