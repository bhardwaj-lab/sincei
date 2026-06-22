from __future__ import annotations

import typer

from ._common_args import (
    BAM_OPTS,
    FILTER_OPTS,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    READ_OPTS,
    DuplicateFilter,
    log_parameters,
    preprocess_args,
)

DESCRIPTION = (
    "Compare read coverages on sampled regions using the Jensen-Shannon distance.\n\n"
    "``scJSD`` samples regions in the genome from BAM files and compares the cumulative read coverages"
    "for each cell on those regions to a synthetic cell with poisson distributed reads using the Jensen-Shannon"
    "distance. Cells with high enrichment of signals show a higher JSD compared to cells whose signal is"
    "homogeneously distributed."
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
    bam_files: list[str] = INPUT_OUTPUT_OPTS["bam_files"],
    barcodes: str = INPUT_OUTPUT_OPTS["barcodes"],
    out_file: str = INPUT_OUTPUT_OPTS["out_file"],
    cell_tag: str = BAM_OPTS["cell_tag"],
    group_tag: str = BAM_OPTS["group_tag"],
    labels: list[str] = BAM_OPTS["labels"],
    smart_labels: bool = BAM_OPTS["smart_labels"],
    blacklist: list[str] = BAM_OPTS["blacklist"],
    chr_to_skip: list[str] = BAM_OPTS["chr_to_skip"],
    bin_size: int = BAM_OPTS["bin_size"],
    duplicate_filter: DuplicateFilter | None = FILTER_OPTS["duplicate_filter"],
    min_mapping_quality: int | None = READ_OPTS["min_mapping_quality"],
    sam_flag_include: int | None = READ_OPTS["sam_flag_include"],
    sam_flag_exclude: int | None = READ_OPTS["sam_flag_exclude"],
    min_fragment_length: int = READ_OPTS["min_fragment_length"],
    max_fragment_length: int = READ_OPTS["max_fragment_length"],
    min_aligned_fraction: float | None = FILTER_OPTS["min_aligned_fraction"],
    number_of_samples: int = typer.Option(
        100000,
        "-n",
        "--numberOfSamples",
        rich_help_panel="Sampling options",
        help="The number of bins that are sampled from the genome, for which the overlapping number of reads is "
        "computed.",
    ),
    skip_zeros: bool = typer.Option(
        False,
        "--skipZeros",
        rich_help_panel="Sampling options",
        help="If set, regions with zero overlapping reads for *all* given BAM files are ignored. This results in a "
        "reduced number of read counts compared to --numberOfSamples.",
    ),
    number_of_processors: str = OTHER_OPTS["number_of_processors"],
    verbose: bool = OTHER_OPTS["verbose"],
    help: bool = OTHER_OPTS["help"],
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
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
