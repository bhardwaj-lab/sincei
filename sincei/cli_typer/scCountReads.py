from __future__ import annotations

import typer

from sincei import _sincei as internal

from . import _parsers as backend
from ._common_args import (
    BAM_OPTS,
    FILTER_OPTS,
    GTF_GFF_OPTS,
    INPUT_OUTPUT_OPTS,
    _IO,
    OTHER_OPTS,
    READ_OPTS,
    Compression,
    DuplicateFilter,
    FilterRNAStrand,
    log_parameters,
    preprocess_args,
)

DESCRIPTION = (
    "Counts reads for each cell barcode on genomic bins or user-defined features.\n\n"
    "``scCountReads`` computes the read coverages per cell barcode for genomic regions in the provided "
    "BAM file(s). The analysis can be performed for the entire genome by running the program in ``bins`` mode. If "
    "you want to count the read coverage for specific regions only, use the ``features`` mode instead. The "
    'standard output of ``scCountReads`` is a ".h5ad" file with counts, along with rowName (features) and colNames '
    "(cell barcodes)."
)


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)

bins_app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help="Count reads in fixed-size bins.",
    context_settings={"help_option_names": []},
)
features_app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help="Count reads on provided features.",
    context_settings={"help_option_names": []},
)

app.add_typer(bins_app, name="bins")
app.add_typer(features_app, name="features")

_COUNTING = "Counting options"

# Tool-specific options reused by both subcommands.
VALUE_TAG = typer.Option(
    None,
    "--valuetag",
    metavar="XX",
    rich_help_panel=_COUNTING,
    help='Instead of counting each read/fragment as "1", add the values from a given BAM tag to the count matrix. '
    "For example, this can be used to count the number of methylated CpG per read.",
)
GENOME_CHUNK_SIZE = typer.Option(
    None,
    "--genomeChunkSize",
    metavar="INT",
    rich_help_panel=_COUNTING,
    help="Manually specify the size of the genome provided to each processor. The default (None) determines this "
    "from the read density of the BAM file.",
)
COMPRESSION = typer.Option(
    Compression.none,
    "--compression",
    metavar="METHOD",
    rich_help_panel=_IO,
    help="HDF5 compression for the output .h5ad datasets.\n\n"
    "[bold yellow]none[/bold yellow]: no compression (default; fastest, larger files).\n\n"
    "[bold yellow]gzip[/bold yellow]: deflate compression (smaller filesize).",
)
COMPRESSION_LEVEL = typer.Option(
    4,
    "--compressionLevel",
    metavar="INT",
    rich_help_panel=_IO,
    help="Compression level (0-9) when --compression is gzip. Ignored for 'none'.",
)


@app.callback(invoke_without_command=True)
def main(ctx: typer.Context, help: bool = OTHER_OPTS["help"]) -> int:
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
        raise typer.Exit()
    return 0


def _count_reads(
    *,
    mode: str,
    bam_files: list[str],
    barcodes: str,
    out_file: str,
    bed: str | None = None,
    labels: list[str] | None,
    smart_labels: bool,
    region: str | None,
    blacklist: list[str] | None,
    chr_to_skip: list[str] | None,
    bin_size: int,
    distance_between_bins: int,
    cell_tag: str,
    group_tag: str | None,
    min_mapping_quality: int | None,
    sam_flag_include: int | None,
    sam_flag_exclude: int | None,
    min_fragment_length: int,
    max_fragment_length: int,
    filter_rna_strand: FilterRNAStrand | None,
    extend_reads: int | None,
    center_reads: bool,
    duplicate_filter: DuplicateFilter | None,
    motif_filter: list[str] | None,
    genome_2bit: str | None,
    gc_content_filter: str | None,
    min_aligned_fraction: float | None,
    value_tag: str | None,
    genome_chunk_size: int | None,
    transcript_id: list[str] | None = None,
    exon_id: list[str] | None = None,
    transcript_id_tag: str | None = None,
    metagene: bool = False,
    compression: Compression = Compression.none,
    compression_level: int = 4,
    number_of_processors: int = 0,
    verbose: bool = False,
) -> int:
    if verbose:
        log_parameters(mode=mode, bam_files=bam_files, out_file=out_file)

    # Options accepted by the CLI for parity with deepTools but not (yet) implemented
    backend.warn_unsupported(
        group_tag=group_tag,
        labels=labels,
        smart_labels=smart_labels,
    )

    min_gc, max_gc = backend.parse_gc_content(gc_content_filter)

    shared = dict(
        barcodes=backend.read_barcodes(barcodes),
        output_path=out_file,
        bc_tag=cell_tag,
        umi_tag=None,
        count_tag=value_tag,
        min_mapq=min_mapping_quality,
        sam_flag_include=sam_flag_include,
        sam_flag_exclude=sam_flag_exclude,
        chr_to_skip=chr_to_skip or [],
        region=region,
        blacklist_path=backend.first_blacklist(blacklist),
        extend_reads=extend_reads,
        center_reads=center_reads,
        dup_method=backend.dup_method(duplicate_filter),
        filter_rna_strand=filter_rna_strand.value if filter_rna_strand else None,
        genome_2bit=genome_2bit,
        motif_filter=backend.parse_motif_filter(motif_filter),
        min_gc=min_gc,
        max_gc=max_gc,
        min_aligned_fraction=min_aligned_fraction,
        min_fragment_length=backend.optional_length(min_fragment_length),
        max_fragment_length=backend.optional_length(max_fragment_length),
        compression=compression.value,
        compression_level=compression_level,
        num_threads=number_of_processors,
    )
    if genome_chunk_size:
        shared["chunk_size"] = genome_chunk_size

    if mode == "bins":
        # stepSize = binSize + gap; a zero gap yields contiguous bins.
        step_size = bin_size + (distance_between_bins or 0)
        internal.count_bins(bam_files, bin_size=bin_size, step_size=step_size, **shared)
    else:
        if not bed:
            raise typer.BadParameter("--bed is required for 'features' mode.")
        # GTF/GFF field selection (ignored for BED inputs). `transcript_id_tag`
        # is None by default so the backend picks the per-format name attribute.
        internal.count_features(
            bam_files,
            bed,
            feature_type=transcript_id,
            exon_type=exon_id,
            name_attr=transcript_id_tag,
            metagene=metagene,
            **shared,
        )

    return 0


@bins_app.callback(invoke_without_command=True)
def bins(
    # Input / Output options
    bam_files: list[str] = INPUT_OUTPUT_OPTS["bam_files"],
    barcodes: str = INPUT_OUTPUT_OPTS["barcodes"],
    out_file: str = INPUT_OUTPUT_OPTS["out_file"],
    region: str | None = INPUT_OUTPUT_OPTS["region"],
    compression: Compression = COMPRESSION,
    compression_level: int = COMPRESSION_LEVEL,
    # BAM options
    labels: list[str] = BAM_OPTS["labels"],
    smart_labels: bool = BAM_OPTS["smart_labels"],
    blacklist: list[str] = BAM_OPTS["blacklist"],
    chr_to_skip: list[str] = BAM_OPTS["chr_to_skip"],
    bin_size: int = BAM_OPTS["bin_size"],
    distance_between_bins: int = BAM_OPTS["distance_between_bins"],
    cell_tag: str = BAM_OPTS["cell_tag"],
    group_tag: str = BAM_OPTS["group_tag"],
    # Read options
    min_mapping_quality: int | None = READ_OPTS["min_mapping_quality"],
    sam_flag_include: int | None = READ_OPTS["sam_flag_include"],
    sam_flag_exclude: int | None = READ_OPTS["sam_flag_exclude"],
    min_fragment_length: int = READ_OPTS["min_fragment_length"],
    max_fragment_length: int = READ_OPTS["max_fragment_length"],
    filter_rna_strand: FilterRNAStrand | None = READ_OPTS["filter_rna_strand"],
    extend_reads: int | None = READ_OPTS["extend_reads"],
    center_reads: bool = READ_OPTS["center_reads"],
    # Filtering options
    duplicate_filter: DuplicateFilter | None = FILTER_OPTS["duplicate_filter"],
    motif_filter: list[str] = FILTER_OPTS["motif_filter"],
    genome_2bit: str | None = FILTER_OPTS["genome_2bit"],
    gc_content_filter: str | None = FILTER_OPTS["gc_content_filter"],
    min_aligned_fraction: float | None = FILTER_OPTS["min_aligned_fraction"],
    # Counting options
    value_tag: str | None = VALUE_TAG,
    genome_chunk_size: int | None = GENOME_CHUNK_SIZE,
    # Other options
    number_of_processors: str = OTHER_OPTS["number_of_processors"],
    verbose: bool = OTHER_OPTS["verbose"],
    help: bool = OTHER_OPTS["help"],
) -> int:
    return _count_reads(
        mode="bins",
        bam_files=bam_files,
        barcodes=barcodes,
        out_file=out_file,
        labels=labels,
        smart_labels=smart_labels,
        region=region,
        blacklist=blacklist,
        chr_to_skip=chr_to_skip,
        bin_size=bin_size,
        distance_between_bins=distance_between_bins,
        cell_tag=cell_tag,
        group_tag=group_tag,
        min_mapping_quality=min_mapping_quality,
        sam_flag_include=sam_flag_include,
        sam_flag_exclude=sam_flag_exclude,
        min_fragment_length=min_fragment_length,
        max_fragment_length=max_fragment_length,
        filter_rna_strand=filter_rna_strand,
        extend_reads=extend_reads,
        center_reads=center_reads,
        duplicate_filter=duplicate_filter,
        motif_filter=motif_filter,
        genome_2bit=genome_2bit,
        gc_content_filter=gc_content_filter,
        min_aligned_fraction=min_aligned_fraction,
        value_tag=value_tag,
        genome_chunk_size=genome_chunk_size,
        compression=compression,
        compression_level=compression_level,
        number_of_processors=number_of_processors,
        verbose=verbose,
    )


@features_app.callback(invoke_without_command=True)
def features(
    # Input / Output options
    bam_files: list[str] = INPUT_OUTPUT_OPTS["bam_files"],
    barcodes: str = INPUT_OUTPUT_OPTS["barcodes"],
    out_file: str = INPUT_OUTPUT_OPTS["out_file"],
    bed: str = typer.Option(
        ...,
        "--bed",
        metavar=".bed/.gtf/.gff",
        rich_help_panel=_IO,
        help="BED/GTF/GFF files to limit the coverage analysis to the regions in them.",
    ),
    region: str | None = INPUT_OUTPUT_OPTS["region"],
    compression: Compression = COMPRESSION,
    compression_level: int = COMPRESSION_LEVEL,
    # BAM options
    labels: list[str] = BAM_OPTS["labels"],
    smart_labels: bool = BAM_OPTS["smart_labels"],
    blacklist: list[str] = BAM_OPTS["blacklist"],
    chr_to_skip: list[str] = BAM_OPTS["chr_to_skip"],
    bin_size: int = BAM_OPTS["bin_size"],
    distance_between_bins: int = BAM_OPTS["distance_between_bins"],
    cell_tag: str = BAM_OPTS["cell_tag"],
    group_tag: str = BAM_OPTS["group_tag"],
    # Read options
    min_mapping_quality: int | None = READ_OPTS["min_mapping_quality"],
    sam_flag_include: int | None = READ_OPTS["sam_flag_include"],
    sam_flag_exclude: int | None = READ_OPTS["sam_flag_exclude"],
    min_fragment_length: int = READ_OPTS["min_fragment_length"],
    max_fragment_length: int = READ_OPTS["max_fragment_length"],
    filter_rna_strand: FilterRNAStrand | None = READ_OPTS["filter_rna_strand"],
    extend_reads: int | None = READ_OPTS["extend_reads"],
    center_reads: bool = READ_OPTS["center_reads"],
    # Filtering options
    duplicate_filter: DuplicateFilter | None = FILTER_OPTS["duplicate_filter"],
    motif_filter: list[str] = FILTER_OPTS["motif_filter"],
    genome_2bit: str | None = FILTER_OPTS["genome_2bit"],
    gc_content_filter: str | None = FILTER_OPTS["gc_content_filter"],
    min_aligned_fraction: float | None = FILTER_OPTS["min_aligned_fraction"],
    # Counting options
    value_tag: str | None = VALUE_TAG,
    genome_chunk_size: int | None = GENOME_CHUNK_SIZE,
    # GTF / GFF options (only affect GTF/GFF inputs, ignored for BED)
    transcript_id: list[str] | None = GTF_GFF_OPTS["transcript_id"],
    exon_id: list[str] | None = GTF_GFF_OPTS["exon_id"],
    transcript_id_tag: str | None = GTF_GFF_OPTS["transcript_id_tag"],
    metagene: bool = typer.Option(
        False,
        "--metagene",
        rich_help_panel="GTF / GFF options",
        help="Group exon intervals by gene (gene_id attribute) and count each read once per gene. "
        "Reads overlapping multiple exons of the same gene are counted once, assigned to the gene "
        "with the greatest total exon overlap.",
    ),
    # Other options
    number_of_processors: str = OTHER_OPTS["number_of_processors"],
    verbose: bool = OTHER_OPTS["verbose"],
    help: bool = OTHER_OPTS["help"],
) -> int:
    return _count_reads(
        mode="features",
        bam_files=bam_files,
        barcodes=barcodes,
        out_file=out_file,
        bed=bed,
        labels=labels,
        smart_labels=smart_labels,
        region=region,
        blacklist=blacklist,
        chr_to_skip=chr_to_skip,
        bin_size=bin_size,
        distance_between_bins=distance_between_bins,
        cell_tag=cell_tag,
        group_tag=group_tag,
        min_mapping_quality=min_mapping_quality,
        sam_flag_include=sam_flag_include,
        sam_flag_exclude=sam_flag_exclude,
        min_fragment_length=min_fragment_length,
        max_fragment_length=max_fragment_length,
        filter_rna_strand=filter_rna_strand,
        extend_reads=extend_reads,
        center_reads=center_reads,
        duplicate_filter=duplicate_filter,
        motif_filter=motif_filter,
        genome_2bit=genome_2bit,
        gc_content_filter=gc_content_filter,
        min_aligned_fraction=min_aligned_fraction,
        value_tag=value_tag,
        genome_chunk_size=genome_chunk_size,
        compression=compression,
        compression_level=compression_level,
        transcript_id=transcript_id,
        exon_id=exon_id,
        transcript_id_tag=transcript_id_tag,
        metagene=metagene,
        number_of_processors=number_of_processors,
        verbose=verbose,
    )


def cli() -> None:
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
