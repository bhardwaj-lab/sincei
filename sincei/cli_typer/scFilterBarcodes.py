from __future__ import annotations

import typer

from sincei import _sincei as internal

from . import _parsers as backend
from ._common_args import (
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
    "Filter cell barcodes from BAM file (for droplet-based single-cell seq).\n\n"
    "``scFilterBarcodes`` identifies barcodes present in a BAM file and produces a list. You can"
    "optionally filter these barcodes by matching them to a whitelist or based on total counts."
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
    bam_file: str = INPUT_OUTPUT_OPTS["bam_file"],
    whitelist: str | None = INPUT_OUTPUT_OPTS["whitelist"],
    out_file: str = INPUT_OUTPUT_OPTS["out_file"],
    region: str | None = INPUT_OUTPUT_OPTS["region"],
    # Barcode options
    min_hamming_dist: int = typer.Option(
        0,
        "-d",
        "--min-hamming-dist",
        metavar="INT",
        rich_help_panel=_BARCODE,
        help="Minimum hamming distance to match the barcode in whitelist. Note that increasing the hamming distance "
        "really slows down the barcode detection process.",
    ),
    min_count: int = typer.Option(
        0,
        "-mc",
        "--min-count",
        metavar="INT",
        rich_help_panel=_BARCODE,
        help="Minimum number of bins with non-zero counts in order to report a barcode. Note that this number ranges "
        "from 0 to genome size / bin size.",
    ),
    rank_plot: str = typer.Option(
        None,
        "-rp",
        "--rank-plot",
        metavar="STR",
        rich_help_panel=_BARCODE,
        help='The output file name to plot the ranked counts per barcode (similar to the "knee plot", but counts are '
        "the number of non-zero bins in this case).",
    ),
    min_mapping_quality: int | None = override(READ_OPTS["min_mapping_quality"], rich_help_panel=_BARCODE),
    # BAM options
    cell_tag: str = BAM_OPTS["cell_tag"],
    group_tag: str = BAM_OPTS["group_tag"],
    labels: list[str] = BAM_OPTS["labels"],
    smart_labels: bool = BAM_OPTS["smart_labels"],
    blacklist: list[str] = BAM_OPTS["blacklist"],
    chr_to_skip: list[str] = BAM_OPTS["chr_to_skip"],
    bin_size: int = override(BAM_OPTS["bin_size"], default=100000),
    distance_between_bins: int | None = override(BAM_OPTS["distance_between_bins"], default=None),
    # Filter options
    duplicate_filter: DuplicateFilter | None = FILTER_OPTS["duplicate_filter"],
    motif_filter: list[str] = FILTER_OPTS["motif_filter"],
    genome_2bit: str | None = FILTER_OPTS["genome_2bit"],
    gc_content_filter: str | None = FILTER_OPTS["gc_content_filter"],
    min_aligned_fraction: float | None = FILTER_OPTS["min_aligned_fraction"],
    # Read options
    min_fragment_length: int = READ_OPTS["min_fragment_length"],
    max_fragment_length: int = READ_OPTS["max_fragment_length"],
    filter_rna_strand: FilterRNAStrand | None = READ_OPTS["filter_rna_strand"],
    extend_reads: int | None = READ_OPTS["extend_reads"],
    center_reads: bool = READ_OPTS["center_reads"],
    # Other options
    number_of_processors: str = OTHER_OPTS["number_of_processors"],
    verbose: bool = OTHER_OPTS["verbose"],
    help: bool = OTHER_OPTS["help"],
) -> int:
    if verbose:
        log_parameters(bam_file=bam_file, whitelist=whitelist, out_file=out_file)

    # Options exposed for parity but not honored by the barcode-detection backend.
    backend.warn_unsupported(
        region=region,
        group_tag=group_tag,
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

    barcode_counts, _selected_barcodes = internal.filter_barcodes(
        bam_file,
        whitelist=backend.read_barcodes(whitelist) if whitelist else None,
        blacklist_file_name=backend.first_blacklist(blacklist),
        cell_tag=cell_tag,
        min_hamming_dist=min_hamming_dist,
        min_count=min_count,
        min_mapping_quality=min_mapping_quality,
        bin_size=bin_size,
        chr_to_skip=chr_to_skip or [],
        num_threads=number_of_processors,
    )

    # Match scFilterBarcodes.py: a TSV with a leading (unnamed) index column plus
    # `barcode`, `count` (number of non-zero bins the barcode was seen in), and
    # `selected` (True when count >= --min-count).  All detected barcodes are
    # listed, sorted by descending count for a knee-plot-friendly ordering.
    barcode_counts.sort(key=lambda bc: (-bc[1], bc[0]))
    with open(out_file, "w") as out:
        out.write("\tbarcode\tcount\tselected\n")
        for i, (barcode, count) in enumerate(barcode_counts):
            selected = count >= min_count
            out.write(f"{i}\t{barcode}\t{count}\t{selected}\n")

    if rank_plot:
        from sincei.plotting._barcode_rank import plot_barcode_rank

        plot_barcode_rank(barcode_counts, rank_plot, min_count=min_count or None)

    return 0


def cli() -> None:
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
