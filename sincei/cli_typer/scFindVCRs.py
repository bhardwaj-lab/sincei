from __future__ import annotations

import typer

from sincei import _sincei as internal

from ._common_args import (
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    log_parameters,
    override,
    preprocess_args,
)

DESCRIPTION = (
    "Call variable chromatin regions (VCRs) from binned chromatin data.\n\n"
    "``scFindVCRs`` calls variable chromatin regions (VCRs) from binned chromatin data. It takes a"
    ".h5ad file containing single-cell genomic signal in bins, and outputs BED files with genome"
    "segmentations for different sensitivities.\n"
    "First, a bin-to-bin correlation matrix is computed for each chromosome.\n"
    "Then, the correlation matrix is turned into a score map by convolving a number of square"
    "Gaussian kernels along its main diagonal. Each kernel has a sigma calculated using a maximum"
    "region size to consider. Each kernel produces a 1-D score for each bin, which are stacked"
    "into a matrix where each row corresponds to a kernel scale and each column to a bin.\n"
    "Finally, the PELT change-point detection algorithm is applied to the score map to identify"
    "regions with distinct correlation patterns. This step depends on a penalty parameter that"
    "controls the number of detected regions."
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
    input: str = INPUT_OUTPUT_OPTS["h5ad_file"],
    region: str | None = INPUT_OUTPUT_OPTS["region"],
    bin_size: int = typer.Option(
        ...,
        "-bs",
        "--binsize",
        metavar="INT",
        rich_help_panel="VCR options",
        help="The size of the bins in the input AnnData object.",
    ),
    max_region_size: int = typer.Option(
        None,
        "-mr",
        "--maxRegionSize",
        metavar="INT",
        rich_help_panel="VCR options",
        help="The maximum region size to be considered, in base pairs. Larger regions may increase compute time. "
        "Defaults to 100 times the bin size.",
    ),
    n_kernels: int = typer.Option(
        20,
        "-nk",
        "--nKernels",
        metavar="INT",
        rich_help_panel="VCR options",
        help="The number of kernels to use for the score map. More kernels generally lead to a better segmentation, "
        "but increase the computational cost.",
    ),
    penalties: list[float] = typer.Option(
        [0.05, 0.1, 0.5],
        "-pen",
        "--penalties",
        metavar="FLOAT",
        rich_help_panel="VCR options",
        help="Penalty value(s) for change-point detection. Higher values result in fewer segments. Multiple values "
        "can be provided (separated by space); each produces a separate set of regions, distinguishable in the output "
        'BED file by filtering on the "score" column.',
    ),
    out_file: str = override(
        INPUT_OUTPUT_OPTS["out_file"],
        metavar=".bed",
        help="Name of the output file (BED format) with the genome segmentation result. The penalty threshold that "
        'defines each segment is saved in the "score" column of the BED file, which can be filtered to obtain '
        "non-overlapping segments.",
    ),
    number_of_processors: str = OTHER_OPTS["number_of_processors"],
    verbose: bool = OTHER_OPTS["verbose"],
    help: bool = OTHER_OPTS["help"],
) -> int:
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
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
