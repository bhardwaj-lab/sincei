from __future__ import annotations

import warnings
from typing import Annotated

import anndata as ad
import mudata as md
import typer

from ._common_args import (
    _IO,
    BAM_OPTS,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    CombineMethod,
    configure_logging,
    log_parameters,
    override,
    preprocess_args,
)
from ._common_args import smart_labels as get_smart_labels
from ._parsers import validate_anndata_list

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
    out_file: Annotated[
        str,
        override(
            INPUT_OUTPUT_OPTS["out_file"],
            metavar=".h5ad/.h5mu",
            help=(
                "The file to write results to. With [bold yellow]multi-sample[/bold "
                "yellow] it is a .h5ad file, which the other sincei tools can use. "
                "With [bold yellow]multi-modal[/bold yellow] it is a .h5mu file, which "
                "only scClusterCells can use, for multi-modal clustering."
            ),
        ),
    ],
    labels: Annotated[
        list[str] | None,
        override(
            BAM_OPTS["labels"],
            rich_help_panel=_IO,
            help=(
                "User defined labels instead of default labels from file names, one "
                "per input file, separated by spaces, e.g. ``--labels sample1 "
                "sample2``. With [bold yellow]multi-modal[/bold yellow] the labels "
                "name the data modalities (e.g. RNA, ATAC)."
            ),
        ),
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
    verbose: Annotated[bool, OTHER_OPTS["verbose"]] = False,
    help: Annotated[bool, OTHER_OPTS["help"]] = False,
) -> int:
    if verbose:
        log_parameters(
            input=input,
            out_file=out_file,
            labels=labels,
            smart_labels=smart_labels,
            method=method,
        )
    else:
        warnings.filterwarnings("ignore")

    if labels and len(labels) != len(input):
        msg = (
            f"got {len(labels)} labels for {len(input)} input files; give one label "
            "per file, or omit --labels to use the file names"
        )
        raise typer.BadParameter(msg, param_hint="--labels")
    names = labels or get_smart_labels(input)

    adatas = validate_anndata_list([ad.read_h5ad(path) for path in input], input)

    if method is CombineMethod.multi_sample:
        for name, adata in zip(names, adatas, strict=True):
            adata.obs_names = [f"{name}_{cell}" for cell in adata.obs_names]
        combined = ad.concat(adatas, merge="first")
        typer.echo(f"Combined cells: {combined.n_obs}")
        typer.echo(f"Combined features: {combined.n_vars}")
        combined.write_h5ad(out_file)
    else:
        mdata = md.MuData(dict(zip(names, adatas, strict=True)))
        typer.echo(f"Combined modalities: {len(mdata.mod)}")
        typer.echo(f"Combined cells: {mdata.n_obs}")
        typer.echo(f"Combined features: {mdata.n_vars}")
        mdata.write_h5mu(out_file)

    return 0


def cli() -> None:
    configure_logging()
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
