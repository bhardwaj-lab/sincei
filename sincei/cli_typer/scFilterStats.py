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

# Column order produced by the Rust `filter_stats` backend (BarcodeStat::to_vec).
_STAT_COLUMNS = [
    "total",
    "filtered",
    "blacklisted",
    "low_mapq",
    "missing_flags",
    "excluded_flags",
    "internal_dupes",
    "external_dupes",
    "singletons",
    "wrong_strand",
    "wrong_motif",
    "wrong_gc",
    "low_aligned_fraction",
]


def _format_row(row) -> list[str]:
    """Format one backend stat row for the TSV.

    Mirrors the original ``sincei``: ``total`` is the absolute number of sampled
    reads, while every other column is reported as a percentage of ``total``
    (``count / total * 100``).  When ``total`` is 0 the percentage is undefined,
    so those cells are left blank (matching the original's NaN -> empty output).
    """
    total = row[0]
    cells = [str(total)]
    for value in row[1:]:
        cells.append("" if total == 0 else str(value / total * 100))
    return cells


DESCRIPTION = (
    "Produce per-cell statistics after filtering reads by user-defined criteria.\n\n"
    "``scFilterStats`` estimates the number of reads that would be filtered given a set of criteria"
    "and prints it to the terminal. Furthermore, it tracks the number of singleton reads."
    "The following metrics will always be tracked regardless of what you specify (the order output also matches this):\n"
    "* Total reads (including unmapped)\n"
    "* Mapped reads\n"
    "* Reads in blacklisted regions (--blackListFileName)\n\n"
    "The following metrics are estimated according to the --binSize and --distanceBetweenBins parameters:\n"
    "* Estimated mapped reads filtered (the total number of mapped reads filtered for any reason)\n"
    "* Alignments with a below threshold MAPQ (--minMappingQuality)\n"
    "* Alignments with at least one missing flag (--samFlagInclude)\n"
    "* Alignments with undesirable flags (--samFlagExclude)\n"
    "* Duplicates determined by sincei (--duplicateFilter)\n"
    "* Duplicates marked externally (e.g., by picard)\n"
    "* Singletons (paired-end reads with only one mate aligning)\n"
    "* Wrong strand (due to --filterRNAstrand)\n\n"
    "The sum of these may be more than the total number of reads. Note that alignments are sampled from"
    "bins of size --binSize spaced --distanceBetweenBins apart."
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
    bin_size: int = override(BAM_OPTS["bin_size"], default=100000),
    distance_between_bins: int | None = override(BAM_OPTS["distance_between_bins"], default=1000000),
    duplicate_filter: DuplicateFilter | None = FILTER_OPTS["duplicate_filter"],
    motif_filter: list[str] = FILTER_OPTS["motif_filter"],
    genome_2bit: str | None = FILTER_OPTS["genome_2bit"],
    gc_content_filter: str | None = FILTER_OPTS["gc_content_filter"],
    min_aligned_fraction: float | None = FILTER_OPTS["min_aligned_fraction"],
    min_mapping_quality: int | None = READ_OPTS["min_mapping_quality"],
    sam_flag_include: int | None = READ_OPTS["sam_flag_include"],
    sam_flag_exclude: int | None = READ_OPTS["sam_flag_exclude"],
    filter_rna_strand: FilterRNAStrand | None = READ_OPTS["filter_rna_strand"],
    labels: list[str] = BAM_OPTS["labels"],
    smart_labels: bool = BAM_OPTS["smart_labels"],
    blacklist: list[str] = BAM_OPTS["blacklist"],
    chr_to_skip: list[str] = BAM_OPTS["chr_to_skip"],
    number_of_processors: str = OTHER_OPTS["number_of_processors"],
    verbose: bool = OTHER_OPTS["verbose"],
    help: bool = OTHER_OPTS["help"],
) -> int:
    if verbose:
        log_parameters(bam_files=bam_files, barcodes=barcodes, out_file=out_file)

    min_gc, max_gc = backend.parse_gc_content(gc_content_filter)
    barcode_list = backend.read_barcodes(barcodes)
    sample_labels = backend.resolve_labels(bam_files, labels, smart_labels)

    kwargs = dict(
        barcodes=barcode_list,
        bc_tag=cell_tag,
        umi_tag=None,
        bin_size=bin_size,
        distance_between_bins=distance_between_bins or 0,
        min_mapq=min_mapping_quality,
        sam_flag_include=sam_flag_include,
        sam_flag_exclude=sam_flag_exclude,
        filter_rna_strand=filter_rna_strand.value if filter_rna_strand else None,
        chr_to_skip=chr_to_skip or [],
        blacklist_path=backend.first_blacklist(blacklist),
        dup_method=backend.dup_method(duplicate_filter),
        genome_2bit=genome_2bit,
        motif_filter=backend.parse_motif_filter(motif_filter),
        min_gc=min_gc,
        max_gc=max_gc,
        min_aligned_fraction=min_aligned_fraction,
        num_threads=number_of_processors,
    )

    # The backend processes one BAM at a time; aggregate the per-sample rows
    # into a single TSV with a leading `sample` column.
    with open(out_file, "w") as out:
        out.write("\t".join(["sample", "barcode", *_STAT_COLUMNS]) + "\n")
        for bam, label in zip(bam_files, sample_labels):
            result_barcodes, stat_rows = internal.filter_stats(bam, **kwargs)
            for barcode, row in zip(result_barcodes, stat_rows):
                out.write("\t".join([label, barcode, *_format_row(row)]) + "\n")

    return 0


def cli() -> None:
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
