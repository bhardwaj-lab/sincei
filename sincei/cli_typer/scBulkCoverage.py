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
    NormalizeUsing,
    OutFileFormat,
    log_parameters,
    override,
    preprocess_args,
)

DESCRIPTION = (
    "Get pseudo-bulk coverage per group using a cell->group mapping (output of scClusterCells).\n\n"
    "``scBulkCoverage`` takes alignments of reads or fragments as input (BAM files), along with cell grouping"
    "information, such as barcode -> batch, or barcode -> cluster, as tsv file, and generates a  coverage track"
    "(bigWig or bedGraph) per group as output. The coverage is calculated as the number of reads/fragments per bin,"
    "where bins are short consecutive counting windows of a defined  size. It is possible to extended/change the"
    "length of the reads to better reflect the actual fragment length. ``scBulkCoverage`` offers normalization per"
    "cluster using different methods."
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
    bam_files: list[str] = INPUT_OUTPUT_OPTS["bam_files"],
    group_info: str = INPUT_OUTPUT_OPTS["group_info"],
    out_prefix: str = INPUT_OUTPUT_OPTS["out_prefix"],
    # Coverage options
    out_file_format: OutFileFormat = typer.Option(
        OutFileFormat.bigwig,
        "-of",
        "--outfile-format",
        metavar="FORMAT",
        rich_help_panel=_COVERAGE,
        help="Output file type.\n\nOne of: [bold yellow]bigwig[/bold yellow], [bold yellow]bedgraph[/bold yellow].",
    ),
    normalize_using: NormalizeUsing = typer.Option(
        NormalizeUsing.cpm,
        "-n",
        "--normalize-using",
        metavar="METHOD",
        rich_help_panel=_COVERAGE,
        help="How to normalize the pseudo-bulk counts.\n\n"
        "[bold yellow]CPM[/bold yellow]: normalize each bin to counts per million mapped reads in that group.\n\n"
        "[bold yellow]RPKM[/bold yellow]: reads per kilobase per million mapped reads in that group.\n\n"
        "[bold yellow]Frequency[/bold yellow]: binarize coverage per bin and normalize to total cells per group.\n\n"
        "[bold yellow]Mean[/bold yellow]: get mean signal per bin across cells in each group.\n\n"
        "[bold yellow]None[/bold yellow]: simply return the sum of coverage per group.",
    ),
    ignore_for_normalization: list[str] = typer.Option(
        None,
        "-ig",
        "--ignore-for-normalization",
        metavar="CHR",
        rich_help_panel=_COVERAGE,
        help="Chromosomes to skip while calculating normalization factors.",
    ),
    normalize_by_reference: str = typer.Option(
        None,
        "-nr",
        "--normalize-by-reference",
        rich_help_panel=_COVERAGE,
        help="NOT IMPLEMENTED: Normalize each group of cells by a reference group (which must be present in the "
        "--group-info file). Note that the --normalize-using method is applied beforehand.",
    ),
    scale_factor: float = typer.Option(
        1.0,
        "--scale-factor",
        rich_help_panel=_COVERAGE,
        help="The computed scaling factor (or 1, if not applicable) will be multiplied by this.",
    ),
    mnase: bool = typer.Option(
        False,
        "--mnase",
        rich_help_panel=_COVERAGE,
        help="Determine nucleosome positions from MNase-seq/CUTnRUN data. Only 3 nucleotides at the center of each "
        "fragment are counted. Only fragment lengths between 130 - 200 bp are considered to avoid dinucleosomes or "
        "other artifacts. *NOTE*: Requires paired-end data. A bin size of 1 is recommended.",
    ),
    offset: list[int] = typer.Option(
        None,
        "--offset",
        rich_help_panel=_COVERAGE,
        help="Uses this offset inside of each read as the signal. This is useful in cases like RiboSeq or GROseq. "
        "Negative values indicate offsets from the end of each read. A value of 1 indicates the first base of the "
        "alignment; -1 is the last base. An offset of 0 is not permitted. If two values are specified, they define a "
        "range of positions.",
    ),
    # BAM options
    cell_tag: str = BAM_OPTS["cell_tag"],
    group_tag: str = BAM_OPTS["group_tag"],
    labels: list[str] = BAM_OPTS["labels"],
    smart_labels: bool = BAM_OPTS["smart_labels"],
    region: str | None = INPUT_OUTPUT_OPTS["region"],
    blacklist: list[str] = BAM_OPTS["blacklist"],
    chr_to_skip: list[str] = BAM_OPTS["chr_to_skip"],
    bin_size: int = override(BAM_OPTS["bin_size"], default=100),
    distance_between_bins: int | None = override(BAM_OPTS["distance_between_bins"], default=None),
    # Filter options
    duplicate_filter: DuplicateFilter | None = FILTER_OPTS["duplicate_filter"],
    motif_filter: list[str] = FILTER_OPTS["motif_filter"],
    genome_2bit: str | None = FILTER_OPTS["genome_2bit"],
    gc_content_filter: str | None = FILTER_OPTS["gc_content_filter"],
    min_aligned_fraction: float | None = FILTER_OPTS["min_aligned_fraction"],
    # Read options
    min_mapping_quality: int | None = READ_OPTS["min_mapping_quality"],
    sam_flag_include: int | None = READ_OPTS["sam_flag_include"],
    sam_flag_exclude: int | None = READ_OPTS["sam_flag_exclude"],
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
    if ctx.invoked_subcommand is None:
        if verbose:
            log_parameters(bam_files=bam_files, group_info=group_info, out_prefix=out_prefix)

        # Options exposed for parity but not honored by the coverage backend.
        backend.warn_unsupported(
            normalize_by_reference=normalize_by_reference,
            group_tag=group_tag,
            region=region,
        )

        min_gc, max_gc = backend.parse_gc_content(gc_content_filter)
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
            umi_tag=None,
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
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
