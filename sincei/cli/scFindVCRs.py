from __future__ import annotations

from typing import Annotated

import anndata as ad
import typer

from sincei.tools.VCRfinder import VCRfinder

from ._common_args import (
    DEFAULT_PROCESSORS,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    configure_logging,
    override,
    preprocess_args,
)
from ._parsers import validate_anndata

DESCRIPTION = (
    "Call variable chromatin regions (VCRs) from binned chromatin data.\n\n"
    "``scFindVCRs`` calls variable chromatin regions (VCRs) from binned chromatin "
    "data. It takes a .h5ad file containing single-cell genomic signal in bins, and "
    "outputs a BED file with the genome segmentation.\n\n"
    "For each chromosome, the counts are binarized: a cell has signal in a bin if its "
    "count is not 0. A bin is kept if at least --minCellPct percent of the cells have "
    "signal in it; the bins at most --slack bins from a kept bin are also kept. The "
    "kept bins are split into signal blocks: runs of bins with signal, with at most "
    "--slack bins without signal between two bins with signal.\n\n"
    "In each block, the co-occurrence of two bins is the number of cells with signal "
    "in both, normalized by the number of cells in each bin (the cosine similarity of "
    "the cells with signal). Each block is then segmented exactly, by dynamic "
    "programming, into runs of bins that co-occur. A penalty for each segment, in "
    "units of bins, controls the number of regions: a higher penalty gives fewer, "
    "longer regions."
)


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)

_VCR = "VCR options"


@app.callback(invoke_without_command=True)
def main(
    input: Annotated[str, INPUT_OUTPUT_OPTS["h5ad_file"]],
    bin_size: Annotated[
        int,
        typer.Option(
            "-bs",
            "--binSize",
            metavar="INT",
            rich_help_panel=_VCR,
            help="The size of the bins in the input AnnData object.",
        ),
    ],
    out_file: Annotated[
        str,
        override(
            INPUT_OUTPUT_OPTS["out_file"],
            metavar=".bed",
            help=(
                "Name of the output file (BED format) with the genome segmentation "
                'result. The "score" column has the cost of each segment (0 for a '
                "region of one bin)."
            ),
        ),
    ],
    region: Annotated[str | None, INPUT_OUTPUT_OPTS["region"]] = None,
    min_cell_pct: Annotated[
        float,
        typer.Option(
            "-mc",
            "--minCellPct",
            metavar="FLOAT",
            rich_help_panel=_VCR,
            help=(
                "Minimum percentage of cells with signal for a bin to be kept. Bins "
                "below it are removed before the signal blocks are made, unless they "
                "are at most --slack bins from a kept bin."
            ),
        ),
    ] = 1.0,
    slack: Annotated[
        int,
        typer.Option(
            "-sl",
            "--slack",
            metavar="INT",
            rich_help_panel=_VCR,
            help=(
                "Number of bins to the left and to the right of each bin that passes "
                "--minCellPct that are also kept, and the maximum number of bins "
                "without signal between two bins with signal in the same block."
            ),
        ),
    ] = 0,
    penalty: Annotated[
        float,
        typer.Option(
            "-pen",
            "--penalty",
            metavar="FLOAT",
            rich_help_panel=_VCR,
            help=(
                "Penalty for each segment, in units of bins. Higher values result in "
                "fewer, longer segments."
            ),
        ),
    ] = 1.0,
    number_of_processors: Annotated[
        int, OTHER_OPTS["number_of_processors"]
    ] = DEFAULT_PROCESSORS,
    verbose: Annotated[bool, OTHER_OPTS["verbose"]] = False,
    help: Annotated[bool, OTHER_OPTS["help"]] = False,
) -> int:
    adata = validate_anndata(ad.read_h5ad(input), input)

    regions = VCRfinder(
        adata=adata,
        binsize=bin_size,
        penalty=penalty,
        X_pct=min_cell_pct,
        k=slack,
        region=region,
        verbose=verbose,
        n_jobs=number_of_processors,
    )

    regions.to_csv(out_file, sep="\t", header=False, index=False)
    return 0


def cli() -> None:
    configure_logging()
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
