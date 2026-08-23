from __future__ import annotations

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
    NormalizeUsing,
    OutFileFormat,
    configure_logging,
    log_parameters,
    preprocess_args,
)

DESCRIPTION = (
    "Get pseudo-bulk coverage per group using a cell->group mapping (output of "
    "scClusterCells).\n\n"
    "``scBulkCoverage`` takes alignments of reads or fragments as input (BAM files), "
    "along with cell grouping information, such as barcode -> batch, or barcode -> "
    "cluster, as tsv file, and generates a  coverage track (bigWig or bedGraph) per "
    "group as output. The coverage is calculated as the number of reads/fragments per "
    "bin, where bins are short consecutive counting windows of a defined  size. It is "
    "possible to extended/change the length of the reads to better reflect the actual "
    "fragment length. ``scBulkCoverage`` offers normalization per cluster using "
    "different methods."
)


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)

_COVERAGE = "Coverage options"


@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    # Input / Output options
    bam_files: Annotated[list[str], INPUT_OUTPUT_OPTS["bam_files"]],
    group_info: Annotated[str, INPUT_OUTPUT_OPTS["group_info"]],
    out_prefix: Annotated[str, INPUT_OUTPUT_OPTS["out_prefix"]],
    # Coverage options
    out_file_format: Annotated[
        OutFileFormat,
        typer.Option(
            "-of",
            "--outFileFormat",
            metavar="FORMAT",
            rich_help_panel=_COVERAGE,
            help=(
                "Output file type.\n\n"
                "One of: [bold yellow]bigwig[/bold yellow], "
                "[bold yellow]bedgraph[/bold yellow]."
            ),
        ),
    ] = OutFileFormat.bigwig,
    normalize_using: Annotated[
        NormalizeUsing,
        typer.Option(
            "-n",
            "--normalizeUsing",
            metavar="METHOD",
            rich_help_panel=_COVERAGE,
            help=(
                "How to normalize the pseudo-bulk counts.\n\n"
                "[bold yellow]CPM[/bold yellow]: normalize each bin to counts per "
                "million mapped reads in that group.\n\n"
                "[bold yellow]RPKM[/bold yellow]: reads per kilobase per million "
                "mapped reads in that group.\n\n"
                "[bold yellow]Frequency[/bold yellow]: binarize coverage per bin and "
                "normalize to total cells per group.\n\n"
                "[bold yellow]Mean[/bold yellow]: get mean signal per bin across cells "
                "in each group.\n\n"
                "[bold yellow]None[/bold yellow]: simply return the sum of coverage "
                "per group."
            ),
        ),
    ] = NormalizeUsing.cpm,
    ignore_for_normalization: Annotated[
        list[str] | None,
        typer.Option(
            "-ig",
            "--ignoreForNormalization",
            metavar="CHR",
            rich_help_panel=_COVERAGE,
            help="Chromosomes to skip while calculating normalization factors.",
        ),
    ] = None,
    scale_factor: Annotated[
        float,
        typer.Option(
            "--scaleFactor",
            rich_help_panel=_COVERAGE,
            help=(
                "The computed scaling factor (or 1, if not applicable) will be "
                "multiplied by this."
            ),
        ),
    ] = 1.0,
    mnase: Annotated[
        bool,
        typer.Option(
            "--mnase",
            rich_help_panel=_COVERAGE,
            help=(
                "Determine nucleosome positions from MNase-seq/CUTnRUN data. Only 3 "
                "nucleotides at the center of each fragment are counted. Only fragment "
                "lengths between 130-200 bp are considered to avoid dinucleosomes or "
                "other artifacts. *NOTE*: Requires paired-end data. A bin size of 1 is "
                "recommended."
            ),
        ),
    ] = False,
    offset: Annotated[
        list[int] | None,
        typer.Option(
            "--offset",
            rich_help_panel=_COVERAGE,
            help=(
                "Uses this offset inside of each read as the signal. This is useful in "
                "cases like RiboSeq or GROseq. Negative values indicate offsets from "
                "the end of each read. A value of 1 indicates the first base of the "
                "alignment; -1 is the last base. An offset of 0 is not permitted. If "
                "two values are specified, they define a range of positions."
            ),
        ),
    ] = None,
    # BAM options
    cell_tag: Annotated[str, BAM_OPTS["cell_tag"]] = "BC",
    umi_tag: Annotated[str, BAM_OPTS["umi_tag"]] = "RX",
    group_tag: Annotated[str | None, BAM_OPTS["group_tag"]] = None,
    labels: Annotated[list[str] | None, BAM_OPTS["labels"]] = None,
    smart_labels: Annotated[bool, BAM_OPTS["smart_labels"]] = False,
    region: Annotated[str | None, INPUT_OUTPUT_OPTS["region"]] = None,
    blacklist: Annotated[list[str] | None, BAM_OPTS["blacklist"]] = None,
    chr_to_skip: Annotated[list[str] | None, BAM_OPTS["chr_to_skip"]] = None,
    bin_size: Annotated[int, BAM_OPTS["bin_size"]] = 100,
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
    min_mapping_quality: Annotated[int | None, READ_OPTS["min_mapping_quality"]] = None,
    sam_flag_include: Annotated[int | None, READ_OPTS["sam_flag_include"]] = None,
    sam_flag_exclude: Annotated[int | None, READ_OPTS["sam_flag_exclude"]] = None,
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
    if ctx.invoked_subcommand is None:
        if verbose:
            log_parameters(
                bam_files=bam_files, group_info=group_info, out_prefix=out_prefix
            )

        backend.require_single_bam_for_group_tag(bam_files, group_tag)

        min_gc, max_gc = backend.parse_gc_content(gc_content_filter)
        if group_tag is not None:
            backend.warn_labels_ignored_under_group_tag(labels, smart_labels)
            # Placeholder only: with --groupTag the group-info file is
            # matched against the BAM's @RG IDs, not against these.
            bam_labels = backend.resolve_labels(bam_files, None, False)
        else:
            bam_labels = backend.resolve_labels(bam_files, labels, smart_labels)
        # stepSize = binSize + gap; a zero gap yields contiguous bins.
        step_size = bin_size + (distance_between_bins or 0)

        output_files = internal.bulk_coverage(
            bam_files,
            bam_labels,
            group_info,
            out_prefix,
            bin_size=bin_size,
            step_size=step_size,
            bc_tag=cell_tag,
            umi_tag=backend.umi_tag_if_used(umi_tag, duplicate_filter),
            group_tag=group_tag,
            region=region,
            min_mapq=min_mapping_quality,
            sam_flag_include=sam_flag_include,
            sam_flag_exclude=sam_flag_exclude,
            chr_to_skip=chr_to_skip or [],
            ignore_for_normalization=ignore_for_normalization or [],
            blacklist_path=backend.first_blacklist(blacklist),
            extend_reads=extend_reads,
            center_reads=center_reads,
            dup_method=backend.dup_method(duplicate_filter),
            genome_2bit=genome_2bit,
            motif_filter=backend.parse_motif_filter(motif_filter),
            min_gc=min_gc,
            max_gc=max_gc,
            min_aligned_fraction=min_aligned_fraction,
            min_fragment_length=backend.optional_length(min_fragment_length),
            max_fragment_length=backend.optional_length(max_fragment_length),
            normalize_using=normalize_using.value,
            scale_factor=scale_factor,
            out_format=out_file_format.value,
            mnase=mnase,
            offset=offset or None,
            filter_rna_strand=filter_rna_strand.value if filter_rna_strand else None,
            num_threads=number_of_processors,
        )
        for path in output_files:
            typer.echo(path)
    return 0


def cli() -> None:
    configure_logging()
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
