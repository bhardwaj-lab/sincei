from __future__ import annotations

from typing import Annotated

import typer

from ._common_args import (
    AVAILABLE_PROCESSORS,
    BAM_OPTS,
    FILTER_OPTS,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    READ_OPTS,
    DuplicateFilter,
    configure_logging,
    log_parameters,
    preprocess_args,
)

DESCRIPTION = (
    "Compare read coverages on sampled regions using the Jensen-Shannon distance.\n\n"
    "``scJSD`` samples regions in the genome from BAM files and compares the "
    "cumulative read coverages for each cell on those regions to a synthetic cell with "
    "poisson distributed reads using the Jensen-Shannon distance. Cells with high "
    "enrichment of signals show a higher JSD compared to cells whose signal is"
    "homogeneously distributed."
)


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)

_SAMPLING = "Sampling options"


@app.callback(invoke_without_command=True)
def main(
    bam_files: Annotated[list[str], INPUT_OUTPUT_OPTS["bam_files"]],
    barcodes: Annotated[str, INPUT_OUTPUT_OPTS["barcodes"]],
    out_file: Annotated[str, INPUT_OUTPUT_OPTS["out_file"]],
    cell_tag: Annotated[str, BAM_OPTS["cell_tag"]] = "BC",
    group_tag: Annotated[str | None, BAM_OPTS["group_tag"]] = None,
    labels: Annotated[list[str] | None, BAM_OPTS["labels"]] = None,
    smart_labels: Annotated[bool, BAM_OPTS["smart_labels"]] = False,
    blacklist: Annotated[list[str] | None, BAM_OPTS["blacklist"]] = None,
    chr_to_skip: Annotated[list[str] | None, BAM_OPTS["chr_to_skip"]] = None,
    bin_size: Annotated[int, BAM_OPTS["bin_size"]] = 10000,
    duplicate_filter: Annotated[
        DuplicateFilter | None, FILTER_OPTS["duplicate_filter"]
    ] = None,
    min_mapping_quality: Annotated[int | None, READ_OPTS["min_mapping_quality"]] = None,
    sam_flag_include: Annotated[int | None, READ_OPTS["sam_flag_include"]] = None,
    sam_flag_exclude: Annotated[int | None, READ_OPTS["sam_flag_exclude"]] = None,
    min_fragment_length: Annotated[int, READ_OPTS["min_fragment_length"]] = 0,
    max_fragment_length: Annotated[int, READ_OPTS["max_fragment_length"]] = 0,
    min_aligned_fraction: Annotated[
        float | None, FILTER_OPTS["min_aligned_fraction"]
    ] = None,
    number_of_samples: Annotated[
        int,
        typer.Option(
            "-n",
            "--numberOfSamples",
            rich_help_panel=_SAMPLING,
            help=(
                "The number of bins that are sampled from the genome, for which the "
                "overlapping number of reads is computed."
            ),
        ),
    ] = 100000,
    skip_zeros: Annotated[
        bool,
        typer.Option(
            "--skipZeros",
            rich_help_panel=_SAMPLING,
            help=(
                "If set, regions with zero overlapping reads for *all* given BAM files "
                "are ignored. This results in a reduced number of read counts compared "
                "to --numberOfSamples."
            ),
        ),
    ] = False,
    number_of_processors: Annotated[
        int, OTHER_OPTS["number_of_processors"]
    ] = AVAILABLE_PROCESSORS,
    verbose: Annotated[bool, OTHER_OPTS["verbose"]] = False,
    help: Annotated[bool, OTHER_OPTS["help"]] = False,
) -> int:
    log_parameters(
        bam_files=bam_files,
        barcodes=barcodes,
        out_file=out_file,
        cell_tag=cell_tag,
        group_tag=group_tag,
        labels=labels,
        smart_labels=smart_labels,
        blacklist=blacklist,
        chr_to_skip=chr_to_skip,
        bin_size=bin_size,
        duplicate_filter=duplicate_filter,
        min_mapping_quality=min_mapping_quality,
        sam_flag_include=sam_flag_include,
        sam_flag_exclude=sam_flag_exclude,
        min_fragment_length=min_fragment_length,
        max_fragment_length=max_fragment_length,
        min_aligned_fraction=min_aligned_fraction,
        number_of_samples=number_of_samples,
        skip_zeros=skip_zeros,
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
