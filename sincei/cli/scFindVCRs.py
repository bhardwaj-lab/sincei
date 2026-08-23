from __future__ import annotations

from typing import Annotated

import typer

from sincei import _sincei as internal

from ._common_args import (
    AVAILABLE_PROCESSORS,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    configure_logging,
    log_parameters,
    override,
    preprocess_args,
)

DESCRIPTION = (
    "Call variable chromatin regions (VCRs) from binned chromatin data.\n\n"
    "``scFindVCRs`` calls variable chromatin regions (VCRs) from binned chromatin "
    "data. It takes a.h5ad file containing single-cell genomic signal in bins, and "
    "outputs BED files with genome segmentations for different sensitivities.\n"
    "First, a bin-to-bin correlation matrix is computed for each chromosome.\n"
    "Then, the correlation matrix is turned into a score map by convolving a number "
    "of square Gaussian kernels along its main diagonal. Each kernel has a sigma "
    "calculated using a maximum region size to consider. Each kernel produces a 1-D "
    "score for each bin, which are stacked into a matrix where each row corresponds to "
    "a kernel scale and each column to a bin.\n"
    "Finally, the PELT change-point detection algorithm is applied to the score map to "
    "identify regions with distinct correlation patterns. This step depends on a "
    "penalty parameter that controls the number of detected regions."
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
            "--binsize",
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
                "result. The penalty threshold that defines each segment is saved in "
                'the "score" column of the BED file, which can be filtered to obtain '
                "non-overlapping segments."
            ),
        ),
    ],
    region: Annotated[str | None, INPUT_OUTPUT_OPTS["region"]] = None,
    max_region_size: Annotated[
        int | None,
        typer.Option(
            "-mr",
            "--maxRegionSize",
            metavar="INT",
            rich_help_panel=_VCR,
            help=(
                "The maximum region size to be considered, in base pairs. Larger "
                "regions may increase compute time. Defaults to 100 times the bin size."
            ),
        ),
    ] = None,
    n_kernels: Annotated[
        int,
        typer.Option(
            "-nk",
            "--nKernels",
            metavar="INT",
            rich_help_panel=_VCR,
            help=(
                "The number of kernels to use for the score map. More kernels "
                "generally lead to a better segmentation, but increase the "
                "computational cost."
            ),
        ),
    ] = 20,
    penalties: Annotated[
        list[float] | None,
        typer.Option(
            "-pen",
            "--penalties",
            metavar="FLOAT",
            rich_help_panel=_VCR,
            # Defaulted in the body rather than here.
            show_default="0.05, 0.1, 0.5",
            help=(
                "Penalty value(s) for change-point detection. Higher values result in "
                "fewer segments. Multiple values can be provided (separated by space); "
                "each produces a separate set of regions, distinguishable in the "
                'output BED file by filtering on the "score" column.'
            ),
        ),
    ] = None,
    number_of_processors: Annotated[
        int, OTHER_OPTS["number_of_processors"]
    ] = AVAILABLE_PROCESSORS,
    verbose: Annotated[bool, OTHER_OPTS["verbose"]] = False,
    help: Annotated[bool, OTHER_OPTS["help"]] = False,
) -> int:
    penalties = penalties or [0.05, 0.1, 0.5]

    if verbose:
        log_parameters(
            input=input,
            region=region,
            bin_size=bin_size,
            max_region_size=max_region_size,
            n_kernels=n_kernels,
            penalties=penalties,
            out_file=out_file,
            number_of_processors=number_of_processors,
            verbose=verbose,
        )

    result_path = internal.find_vcrs(
        input,
        bin_size,
        out_file,
        max_region_size=max_region_size,
        n_kernels=n_kernels,
        penalties=penalties,
        region=region,
        num_threads=number_of_processors,
    )
    typer.echo(result_path)
    return 0


def cli() -> None:
    configure_logging()
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
