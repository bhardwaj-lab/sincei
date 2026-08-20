from __future__ import annotations

import math
from pathlib import Path
from typing import Annotated

import typer

from sincei import _sincei as internal

from . import _parsers as backend
from ._common_args import (
    AVAILABLE_PROCESSORS,
    BAM_OPTS,
    FILTER_OPTS,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    READ_OPTS,
    DuplicateFilter,
    FilterRNAStrand,
    log_parameters,
    override,
    preprocess_args,
)

DESCRIPTION = (
    "Filter cell barcodes from a BAM file (for droplet-based single-cell seq).\n\n"
    "``scFilterBarcodes`` identifies barcodes present in a BAM file and produces a "
    "list. You can optionally filter these barcodes by matching them to a whitelist or "
    "based on total counts. This tool expects single experiment BAM files, not merged "
    "files."
)


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)

_BARCODE = "Barcode options "


@app.callback(invoke_without_command=True)
def main(
    # Input / Output options
    bam_file: Annotated[str, INPUT_OUTPUT_OPTS["bam_file"]],
    out_file: Annotated[str, INPUT_OUTPUT_OPTS["out_file"]],
    whitelist: Annotated[str | None, INPUT_OUTPUT_OPTS["whitelist"]] = None,
    region: Annotated[str | None, INPUT_OUTPUT_OPTS["region"]] = None,
    # Barcode options
    min_hamming_dist: Annotated[
        int,
        typer.Option(
            "-d",
            "--minHammingDist",
            metavar="INT",
            rich_help_panel=_BARCODE,
            help=(
                "Minimum hamming distance to match the barcode in whitelist. Note that "
                "increasing the hamming distance really slows down the barcode "
                "detection process."
            ),
        ),
    ] = 0,
    min_count: Annotated[
        int,
        typer.Option(
            "-mc",
            "--minCount",
            metavar="INT",
            rich_help_panel=_BARCODE,
            help=(
                "Minimum number of bins with non-zero counts in order to report a "
                "barcode. Note that this number ranges from 0 to genome size / bin "
                "size."
            ),
        ),
    ] = 0,
    rank_plot: Annotated[
        str | None,
        typer.Option(
            "-rp",
            "--rankPlot",
            metavar="STR",
            rich_help_panel=_BARCODE,
            help=(
                "The output file name to plot the ranked counts per barcode (similar "
                'to the "knee plot", but counts are the number of non-zero bins in '
                "this case)."
            ),
        ),
    ] = None,
    min_mapping_quality: Annotated[
        int | None,
        override(READ_OPTS["min_mapping_quality"], rich_help_panel=_BARCODE),
    ] = None,
    # BAM options
    cell_tag: Annotated[str, BAM_OPTS["cell_tag"]] = "BC",
    group_tag: Annotated[str | None, BAM_OPTS["group_tag"]] = None,
    labels: Annotated[list[str] | None, BAM_OPTS["labels"]] = None,
    smart_labels: Annotated[bool, BAM_OPTS["smart_labels"]] = False,
    blacklist: Annotated[list[str] | None, BAM_OPTS["blacklist"]] = None,
    chr_to_skip: Annotated[list[str] | None, BAM_OPTS["chr_to_skip"]] = None,
    bin_size: Annotated[int, BAM_OPTS["bin_size"]] = 100000,
    distance_between_bins: Annotated[
        int | None, BAM_OPTS["distance_between_bins"]
    ] = None,
    # Filter options
    duplicate_filter: Annotated[
        DuplicateFilter | None, FILTER_OPTS["duplicate_filter"]
    ] = None,
    motif_filter: Annotated[list[str] | None, FILTER_OPTS["motif_filter"]] = None,
    genome_2bit: Annotated[str | None, FILTER_OPTS["genome_2bit"]] = None,
    gc_content_filter: Annotated[str | None, FILTER_OPTS["gc_content_filter"]] = None,
    min_aligned_fraction: Annotated[
        float | None, FILTER_OPTS["min_aligned_fraction"]
    ] = None,
    # Read options
    min_fragment_length: Annotated[int, READ_OPTS["min_fragment_length"]] = 0,
    max_fragment_length: Annotated[int, READ_OPTS["max_fragment_length"]] = 0,
    filter_rna_strand: Annotated[
        FilterRNAStrand | None, READ_OPTS["filter_rna_strand"]
    ] = None,
    extend_reads: Annotated[int | None, READ_OPTS["extend_reads"]] = None,
    center_reads: Annotated[bool, READ_OPTS["center_reads"]] = False,
    # Other options
    number_of_processors: Annotated[
        int, OTHER_OPTS["number_of_processors"]
    ] = AVAILABLE_PROCESSORS,
    verbose: Annotated[bool, OTHER_OPTS["verbose"]] = False,
    help: Annotated[bool, OTHER_OPTS["help"]] = False,
) -> int:
    if verbose:
        log_parameters(bam_file=bam_file, whitelist=whitelist, out_file=out_file)

    # Options exposed for parity but not honored by the barcode-detection backend.
    backend.warn_unsupported(
        region=region,
        labels=labels,
        smart_labels=smart_labels,
        distance_between_bins=distance_between_bins,
        duplicate_filter=duplicate_filter,
        motif_filter=motif_filter,
        genome_2bit=genome_2bit,
        gc_content_filter=gc_content_filter,
        min_aligned_fraction=min_aligned_fraction,
        min_fragment_length=backend.optional_length(min_fragment_length),
        max_fragment_length=backend.optional_length(max_fragment_length),
        filter_rna_strand=filter_rna_strand,
        extend_reads=extend_reads,
        center_reads=center_reads,
    )

    barcode_counts = internal.filter_barcodes(
        bam_file,
        whitelist=backend.read_barcodes(whitelist) if whitelist else None,
        blacklist_file_name=backend.first_blacklist(blacklist),
        cell_tag=cell_tag,
        group_tag=group_tag,
        min_hamming_dist=min_hamming_dist,
        min_mapping_quality=min_mapping_quality,
        bin_size=bin_size,
        chr_to_skip=chr_to_skip or [],
        num_threads=number_of_processors,
    )

    # A TSV of every detected barcode, sorted by descending count for a
    # knee-plot-friendly ordering: `count` is the number of non-zero bins the
    # barcode was seen in, `selected` is True when count >= min_count, and
    # `count_log10` and `count_rank` are the knee-plot axes.
    #
    # Two barcodes with the same count share the same rank, and the next rank is
    # incremented by the number of barcodes with that count.
    rank_of_count: dict[int, int] = {}
    for position, (_, count) in enumerate(barcode_counts):
        rank_of_count.setdefault(count, position + 1)

    with Path(out_file).open("w") as out:
        out.write("barcode\tcount\tselected\tcount_log10\tcount_rank\n")
        for barcode, count in barcode_counts:
            selected = count >= min_count
            # The backend only reports barcodes it saw, so count >= 1 and
            # log10 is always defined.
            out.write(
                f"{barcode}\t{count}\t{selected}\t{math.log10(count)}\t"
                f"{rank_of_count[count]}\n"
            )

    if rank_plot:
        # Imported lazily: it pulls in matplotlib, which would otherwise be paid
        # for on every invocation of this command.
        from sincei.plotting._barcode_rank import plot_barcode_rank  # noqa: PLC0415

        plot_barcode_rank(barcode_counts, rank_plot, min_count=min_count or None)

    return 0


def cli() -> None:
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
