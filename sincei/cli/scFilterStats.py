from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Annotated

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
    preprocess_args,
)

if TYPE_CHECKING:
    from collections.abc import Sequence

# Column order produced by the Rust `filter_stats` backend (BarcodeStat::to_vec).
_STAT_COLUMNS = [
    "Total_sampled",
    "Filtered",
    "Blacklisted",
    "Low_MAPQ",
    "Missing_Flags",
    "Excluded_Flags",
    "Internal_Duplicates",
    "Marked_Duplicates",
    "Singletons",
    "Wrong_strand",
    "Wrong_motif",
    "Unwanted_GC_content",
    "Low_aligned_fraction",
]


def _format_row(row: Sequence[float]) -> list[str]:
    """Format one backend stat row for the TSV.

    Mirrors the original ``sincei``: ``total`` is the absolute number of sampled
    reads, while every other column is reported as a percentage of ``total``
    (``count / total * 100``).  When ``total`` is 0 the percentage is undefined,
    so those cells are left blank (matching the original's NaN -> empty output).
    """
    total = row[0]
    return [
        str(total),
        *("" if total == 0 else str(value / total * 100) for value in row[1:]),
    ]


DESCRIPTION = (
    "Produce per-cell statistics after filtering reads by user-defined criteria.\n\n"
    "``scFilterStats`` estimates the number of reads that would be filtered given a "
    "set of criteria and prints it to the terminal. Furthermore, it tracks the number "
    "of singleton reads. The following metrics will always be tracked regardless of "
    "what you specify (the order output also matches this):\n"
    "* Total sampled reads (including unmapped)\n"
    "* Mapped reads\n"
    "* Reads in blacklisted regions (--blackListFileName)\n\n"
    "The following metrics are estimated according to the --binSize and "
    "--distanceBetweenBins parameters:\n"
    "* Estimated mapped reads filtered (the total number of mapped reads filtered for "
    "any reason)\n"
    "* Alignments with a below threshold MAPQ (--minMappingQuality)\n"
    "* Alignments with at least one missing flag (--samFlagInclude)\n"
    "* Alignments with undesirable flags (--samFlagExclude)\n"
    "* Duplicates determined by sincei (--duplicateFilter)\n"
    "* Duplicates marked externally (e.g., by picard)\n"
    "* Singletons (paired-end reads with only one mate aligning)\n"
    "* Wrong strand (due to --filterRNAstrand)\n\n"
    "The sum of these may be more than the total number of reads. Note that alignments "
    "are sampled from bins of size --binSize spaced --distanceBetweenBins apart."
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
    bam_files: Annotated[list[str], INPUT_OUTPUT_OPTS["bam_files"]],
    barcodes: Annotated[str, INPUT_OUTPUT_OPTS["barcodes"]],
    out_file: Annotated[str, INPUT_OUTPUT_OPTS["out_file"]],
    cell_tag: Annotated[str, BAM_OPTS["cell_tag"]] = "BC",
    umi_tag: Annotated[str, BAM_OPTS["umi_tag"]] = "RX",
    bin_size: Annotated[int, BAM_OPTS["bin_size"]] = 100_000,
    distance_between_bins: Annotated[
        int | None, BAM_OPTS["distance_between_bins"]
    ] = 1_000_000,
    duplicate_filter: Annotated[
        DuplicateFilter | None, FILTER_OPTS["duplicate_filter"]
    ] = None,
    motif_filter: Annotated[list[str] | None, FILTER_OPTS["motif_filter"]] = None,
    genome_2bit: Annotated[str | None, FILTER_OPTS["genome_2bit"]] = None,
    gc_content_filter: Annotated[str | None, FILTER_OPTS["gc_content_filter"]] = None,
    min_aligned_fraction: Annotated[
        float | None, FILTER_OPTS["min_aligned_fraction"]
    ] = None,
    min_mapping_quality: Annotated[int | None, READ_OPTS["min_mapping_quality"]] = None,
    sam_flag_include: Annotated[int | None, READ_OPTS["sam_flag_include"]] = None,
    sam_flag_exclude: Annotated[int | None, READ_OPTS["sam_flag_exclude"]] = None,
    filter_rna_strand: Annotated[
        FilterRNAStrand | None, READ_OPTS["filter_rna_strand"]
    ] = None,
    group_tag: Annotated[str | None, BAM_OPTS["group_tag"]] = None,
    labels: Annotated[list[str] | None, BAM_OPTS["labels"]] = None,
    smart_labels: Annotated[bool, BAM_OPTS["smart_labels"]] = False,
    blacklist: Annotated[list[str] | None, BAM_OPTS["blacklist"]] = None,
    chr_to_skip: Annotated[list[str] | None, BAM_OPTS["chr_to_skip"]] = None,
    number_of_processors: Annotated[
        int, OTHER_OPTS["number_of_processors"]
    ] = AVAILABLE_PROCESSORS,
    verbose: Annotated[bool, OTHER_OPTS["verbose"]] = False,
    help: Annotated[bool, OTHER_OPTS["help"]] = False,
) -> int:
    if verbose:
        log_parameters(bam_files=bam_files, barcodes=barcodes, out_file=out_file)

    min_gc, max_gc = backend.parse_gc_content(gc_content_filter)
    barcode_list = backend.read_barcodes(barcodes)
    # A merged BAM carries its samples in a read tag rather than in separate
    # files, so the barcode alone no longer identifies a cell.
    if group_tag is not None and len(bam_files) > 1:
        msg = (
            f"--groupTag expects a single merged BAM, but {len(bam_files)} were "
            f"given. Merge them first, e.g. `samtools merge -r`."
        )
        raise typer.BadParameter(msg)

    sample_labels = backend.resolve_labels(bam_files, labels, smart_labels)

    kwargs = {
        "barcodes": barcode_list,
        "bc_tag": cell_tag,
        "group_tag": group_tag,
        "umi_tag": backend.umi_tag_if_used(umi_tag, duplicate_filter),
        "bin_size": bin_size,
        "distance_between_bins": distance_between_bins or 0,
        "min_mapq": min_mapping_quality,
        "sam_flag_include": sam_flag_include,
        "sam_flag_exclude": sam_flag_exclude,
        "filter_rna_strand": filter_rna_strand.value if filter_rna_strand else None,
        "chr_to_skip": chr_to_skip or [],
        "blacklist_path": backend.first_blacklist(blacklist),
        "dup_method": backend.dup_method(duplicate_filter),
        "genome_2bit": genome_2bit,
        "motif_filter": backend.parse_motif_filter(motif_filter),
        "min_gc": min_gc,
        "max_gc": max_gc,
        "min_aligned_fraction": min_aligned_fraction,
        "num_threads": number_of_processors,
    }

    # The backend processes one BAM at a time; aggregate the per-sample rows
    # into a single TSV keyed by `Cell_ID` (`{sample}::{barcode}`, the same
    # identifier the counting tools use for AnnData's obs_names).
    with Path(out_file).open("w") as out:
        out.write("\t".join(["Cell_ID", *_STAT_COLUMNS]) + "\n")
        for bam, label in zip(bam_files, sample_labels, strict=True):
            row_labels, stat_rows = internal.filter_stats(bam, **kwargs)
            # Grouping makes the backend return `group::barcode` already, since
            # the sample comes from the reads rather than from the file.
            cell_ids = (
                row_labels
                if group_tag is not None
                else [f"{label}::{barcode}" for barcode in row_labels]
            )
            out.writelines(
                "\t".join([cell_id, *_format_row(row)]) + "\n"
                for cell_id, row in zip(cell_ids, stat_rows, strict=True)
            )

    return 0


def cli() -> None:
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
