from __future__ import annotations

from typing import Annotated

import typer

from sincei import _sincei as internal

from . import _parsers as backend
from ._common_args import (
    _IO,
    AVAILABLE_PROCESSORS,
    BAM_OPTS,
    FILTER_OPTS,
    GTF_GFF_OPTS,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    READ_OPTS,
    Compression,
    DuplicateFilter,
    FilterRNAStrand,
    configure_logging,
    log_parameters,
    preprocess_args,
)

DESCRIPTION = (
    "Counts reads for each cell barcode on genomic bins or user-defined features.\n\n"
    "``scCountReads`` computes the read coverages per cell barcode for genomic regions "
    "in the provided BAM file(s). The analysis can be performed for the entire genome "
    "by running the program in ``bins`` mode. If you want to count the read coverage "
    "for specific regions only, use the ``features`` mode instead. The standard output "
    'of ``scCountReads`` is a ".h5ad" file with counts, along with rowName (features) '
    "and colNames (cell barcodes)."
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
    help=(
        "Count reads in fixed-size bins. A read is counted in every bin it "
        "overlaps, so a read on a bin boundary adds to both bins."
    ),
    context_settings={"help_option_names": []},
)
features_app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=(
        "Count reads on provided features. A read is counted once for every "
        "feature it overlaps, so overlapping features each get their own count. "
        "Use --metagene to count each read only once."
    ),
    context_settings={"help_option_names": []},
)

app.add_typer(bins_app, name="bins")
app.add_typer(features_app, name="features")

_COUNTING = "Counting options"

# Tool-specific options reused by both subcommands.
VALUE_TAG = typer.Option(
    "--valuetag",
    metavar="XX",
    rich_help_panel=_COUNTING,
    help=(
        'Instead of counting each read/fragment as "1", add the values from a given '
        "BAM tag to the count matrix. For example, this can be used to count the "
        "number of methylated CpG per read. A read that does not carry the tag is not "
        "counted at all. Negative values contribute their magnitude, so a tag of -3 "
        "adds 3."
    ),
)
GENOME_CHUNK_SIZE = typer.Option(
    "--genomeChunkSize",
    metavar="INT",
    rich_help_panel=_COUNTING,
    help=(
        "Manually specify the size of the genome provided to each processor. "
        "(Default: 1Mb)"
    ),
)
COMPRESSION = typer.Option(
    "--compression",
    metavar="METHOD",
    rich_help_panel=_IO,
    help=(
        "HDF5 compression for the output .h5ad datasets.\n\n"
        "[bold yellow]none[/bold yellow]: no compression (default; fastest, larger "
        "files).\n\n"
        "[bold yellow]gzip[/bold yellow]: deflate compression (smaller filesize)."
    ),
)
COMPRESSION_LEVEL = typer.Option(
    "--compressionLevel",
    metavar="INT",
    rich_help_panel=_IO,
    help="Compression level (0-9) when --compression is gzip. Ignored for 'none'.",
)
BED = typer.Option(
    "--bed",
    metavar=".bed/.gtf/.gff",
    rich_help_panel=_IO,
    help="BED/GTF/GFF files to limit the coverage analysis to the regions in them.",
)


@app.callback(invoke_without_command=True)
def main(ctx: typer.Context, help: Annotated[bool, OTHER_OPTS["help"]] = False) -> int:
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
    # Bins mode only; `features` tiles nothing, so it passes neither and the
    # defaults stand. A zero bin size is not a usable tiling, and the backend
    # rejects it (`bin_matrix.rs`), so a bins caller that forgot it fails loudly.
    bin_size: int = 0,
    distance_between_bins: int = 0,
    cell_tag: str,
    umi_tag: str,
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
    feature_id: list[str] | None = None,
    exon_id: list[str] | None = None,
    feature_id_tag: str | None = None,
    metagene: bool = False,
    compression: Compression = Compression.none,
    compression_level: int = 4,
    number_of_processors: int = 0,
    verbose: bool = False,
) -> int:
    if verbose:
        log_parameters(mode=mode, bam_files=bam_files, out_file=out_file)

    backend.require_single_bam_for_group_tag(bam_files, group_tag)

    # Row labels: the file stem by default (which is also what --smartLabels
    # asks for), or --labels when given. Under --groupTag they come from the
    # BAM's @RG IDs instead, so anything passed here would be ignored.
    if group_tag is not None:
        backend.warn_labels_ignored_under_group_tag(labels, smart_labels)
        sample_labels: list[str] = []
    else:
        sample_labels = backend.resolve_labels(bam_files, labels, smart_labels)

    min_gc, max_gc = backend.parse_gc_content(gc_content_filter)

    shared = {
        "barcodes": backend.read_barcodes(barcodes),
        "labels": sample_labels,
        "output_path": out_file,
        "bc_tag": cell_tag,
        "umi_tag": backend.umi_tag_if_used(umi_tag, duplicate_filter),
        "count_tag": value_tag,
        "group_tag": group_tag,
        "min_mapq": min_mapping_quality,
        "sam_flag_include": sam_flag_include,
        "sam_flag_exclude": sam_flag_exclude,
        "chr_to_skip": chr_to_skip or [],
        "region": region,
        "blacklist_path": backend.first_blacklist(blacklist),
        "extend_reads": extend_reads,
        "center_reads": center_reads,
        "dup_method": backend.dup_method(duplicate_filter),
        "filter_rna_strand": filter_rna_strand.value if filter_rna_strand else None,
        "genome_2bit": genome_2bit,
        "motif_filter": backend.parse_motif_filter(motif_filter),
        "min_gc": min_gc,
        "max_gc": max_gc,
        "min_aligned_fraction": min_aligned_fraction,
        "min_fragment_length": backend.optional_length(min_fragment_length),
        "max_fragment_length": backend.optional_length(max_fragment_length),
        "compression": compression.value,
        "compression_level": compression_level,
        "num_threads": number_of_processors,
    }
    if genome_chunk_size:
        shared["chunk_size"] = genome_chunk_size

    if mode == "bins":
        # stepSize = binSize + gap; a zero gap yields contiguous bins.
        step_size = bin_size + distance_between_bins
        internal.count_bins(bam_files, bin_size=bin_size, step_size=step_size, **shared)
    else:
        if not bed:
            msg = "--bed is required for 'features' mode."
            raise typer.BadParameter(msg)
        # GTF/GFF field selection (ignored for BED inputs). `feature_id_tag`
        # is None by default so the backend picks the per-format name attribute.
        internal.count_features(
            bam_files,
            bed,
            feature_type=feature_id,
            exon_type=exon_id,
            name_attr=feature_id_tag,
            metagene=metagene,
            **shared,
        )

    return 0


@bins_app.callback(invoke_without_command=True)
def bins(
    # Input / Output options
    bam_files: Annotated[list[str], INPUT_OUTPUT_OPTS["bam_files"]],
    barcodes: Annotated[str, INPUT_OUTPUT_OPTS["barcodes"]],
    out_file: Annotated[str, INPUT_OUTPUT_OPTS["out_file"]],
    region: Annotated[str | None, INPUT_OUTPUT_OPTS["region"]] = None,
    compression: Annotated[Compression, COMPRESSION] = Compression.none,
    compression_level: Annotated[int, COMPRESSION_LEVEL] = 4,
    # BAM options
    labels: Annotated[list[str] | None, BAM_OPTS["labels"]] = None,
    smart_labels: Annotated[bool, BAM_OPTS["smart_labels"]] = False,
    blacklist: Annotated[list[str] | None, BAM_OPTS["blacklist"]] = None,
    chr_to_skip: Annotated[list[str] | None, BAM_OPTS["chr_to_skip"]] = None,
    bin_size: Annotated[int, BAM_OPTS["bin_size"]] = 10000,
    distance_between_bins: Annotated[int, BAM_OPTS["distance_between_bins"]] = 0,
    cell_tag: Annotated[str, BAM_OPTS["cell_tag"]] = "BC",
    umi_tag: Annotated[str, BAM_OPTS["umi_tag"]] = "RX",
    group_tag: Annotated[str | None, BAM_OPTS["group_tag"]] = None,
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
    # Filtering options
    duplicate_filter: Annotated[
        DuplicateFilter | None, FILTER_OPTS["duplicate_filter"]
    ] = None,
    motif_filter: Annotated[list[str] | None, FILTER_OPTS["motif_filter"]] = None,
    genome_2bit: Annotated[str | None, FILTER_OPTS["genome_2bit"]] = None,
    gc_content_filter: Annotated[str | None, FILTER_OPTS["gc_content_filter"]] = None,
    min_aligned_fraction: Annotated[
        float | None, FILTER_OPTS["min_aligned_fraction"]
    ] = None,
    # Counting options
    value_tag: Annotated[str | None, VALUE_TAG] = None,
    genome_chunk_size: Annotated[int | None, GENOME_CHUNK_SIZE] = None,
    # Other options
    number_of_processors: Annotated[
        int, OTHER_OPTS["number_of_processors"]
    ] = AVAILABLE_PROCESSORS,
    verbose: Annotated[bool, OTHER_OPTS["verbose"]] = False,
    help: Annotated[bool, OTHER_OPTS["help"]] = False,
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
        umi_tag=umi_tag,
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
    bam_files: Annotated[list[str], INPUT_OUTPUT_OPTS["bam_files"]],
    barcodes: Annotated[str, INPUT_OUTPUT_OPTS["barcodes"]],
    out_file: Annotated[str, INPUT_OUTPUT_OPTS["out_file"]],
    bed: Annotated[str, BED],
    region: Annotated[str | None, INPUT_OUTPUT_OPTS["region"]] = None,
    compression: Annotated[Compression, COMPRESSION] = Compression.none,
    compression_level: Annotated[int, COMPRESSION_LEVEL] = 4,
    # BAM options
    labels: Annotated[list[str] | None, BAM_OPTS["labels"]] = None,
    smart_labels: Annotated[bool, BAM_OPTS["smart_labels"]] = False,
    blacklist: Annotated[list[str] | None, BAM_OPTS["blacklist"]] = None,
    chr_to_skip: Annotated[list[str] | None, BAM_OPTS["chr_to_skip"]] = None,
    cell_tag: Annotated[str, BAM_OPTS["cell_tag"]] = "BC",
    umi_tag: Annotated[str, BAM_OPTS["umi_tag"]] = "RX",
    group_tag: Annotated[str | None, BAM_OPTS["group_tag"]] = None,
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
    # Filtering options
    duplicate_filter: Annotated[
        DuplicateFilter | None, FILTER_OPTS["duplicate_filter"]
    ] = None,
    motif_filter: Annotated[list[str] | None, FILTER_OPTS["motif_filter"]] = None,
    genome_2bit: Annotated[str | None, FILTER_OPTS["genome_2bit"]] = None,
    gc_content_filter: Annotated[str | None, FILTER_OPTS["gc_content_filter"]] = None,
    min_aligned_fraction: Annotated[
        float | None, FILTER_OPTS["min_aligned_fraction"]
    ] = None,
    # Counting options
    value_tag: Annotated[str | None, VALUE_TAG] = None,
    genome_chunk_size: Annotated[int | None, GENOME_CHUNK_SIZE] = None,
    # GTF / GFF options (only affect GTF/GFF inputs, ignored for BED)
    feature_id: Annotated[list[str] | None, GTF_GFF_OPTS["feature_id"]] = None,
    exon_id: Annotated[list[str] | None, GTF_GFF_OPTS["exon_id"]] = None,
    feature_id_tag: Annotated[str | None, GTF_GFF_OPTS["feature_id_tag"]] = None,
    metagene: Annotated[bool, GTF_GFF_OPTS["metagene"]] = False,
    # Other options
    number_of_processors: Annotated[
        int, OTHER_OPTS["number_of_processors"]
    ] = AVAILABLE_PROCESSORS,
    verbose: Annotated[bool, OTHER_OPTS["verbose"]] = False,
    help: Annotated[bool, OTHER_OPTS["help"]] = False,
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
        cell_tag=cell_tag,
        umi_tag=umi_tag,
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
        feature_id=feature_id,
        exon_id=exon_id,
        feature_id_tag=feature_id_tag,
        metagene=metagene,
        number_of_processors=number_of_processors,
        verbose=verbose,
    )


def cli() -> None:
    configure_logging()
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
