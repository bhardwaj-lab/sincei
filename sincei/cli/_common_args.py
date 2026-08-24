"""Shared CLI building blocks for the sincei typer apps.

This module is the single source of truth for:

* the enums used to validate choice-style options;
* the dictionaries of reusable ``typer.Option`` definitions;
* the option parsers/callbacks (``normalize_*``, ``version_callback``,
  ``help_callback``) and the ``override`` helper used to build option variants.

The cross-command plumbing (``log_parameters``, ``preprocess_args``) lives in
``_parsers`` and is re-exported here so existing imports keep working.

Command modules import the option dictionaries and attach the entries with
``Annotated[<type>, <option>]``, so the flags, help text and metavars only have to
be maintained in one place. Options whose *metadata* differs between commands are
reused with ``override(...)``; options whose *default* differs simply get a
different default in the function signature.

The options here are written in typer's ``Annotated`` style, which means they carry
no default value: ``typer.Option`` receives only the flags and the display metadata,
and every default lives in the command signature. That keeps each parameter's
annotation honest about the value the function receives, which is what lets ``ty``
check these modules.
"""

from __future__ import annotations

import copy
import os
from enum import Enum
from pathlib import Path

import click
import typer

from sincei import _sincei as internal

# Re-exported so ``from ._common_args import log_parameters, preprocess_args``
# keeps working for the command modules; the definitions are in ``_parsers``.
from ._parsers import configure_logging, log_parameters, preprocess_args  # noqa: F401


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def version_string() -> str:
    return internal.version()


def normalize_region(
    ctx: click.Context, param: click.Parameter, string: str | None
) -> str | None:
    if string is None:
        return None

    # Remove whitespaces using the split/join trick.
    region = "".join(string.split())
    if region == "":
        return None

    # Drop undesired characters that may be present and replace "-" with ":".
    region = region.translate({
        ord(character): None for character in ",;|!{}()"
    }).replace("-", ":")

    if len(region) == 0:
        msg = f"{string} is not a valid region"
        raise click.BadParameter(msg, param=param, ctx=ctx)

    return region


def available_processors() -> int:
    """Usable processor count (``os.sched_getaffinity`` where it exists)."""
    # Use os.sched_getaffinity if available (Linux), fall back to os.cpu_count()
    # (MacOS/Windows)
    try:
        return len(os.sched_getaffinity(0))  # ty: ignore[unresolved-attribute]
    except (AttributeError, OSError):
        return os.cpu_count() or 1


AVAILABLE_PROCESSORS = available_processors()


def normalize_processors(value: str | int) -> int:
    """Parse ``--numberOfProcessors`` into a concrete processor count.

    Used as the option's ``parser`` rather than its ``callback`` so that click
    does not coerce ``"max"`` / ``"max/2"`` to ``int`` before we see it, which is
    what a plain ``int`` annotation would otherwise make it do.  Typer also runs
    the parser over the default, hence the ``int`` half of the input type.
    """
    if value == "max/2":
        return int(AVAILABLE_PROCESSORS * 0.5)
    if value == "max":
        return AVAILABLE_PROCESSORS

    try:
        number_of_processors = int(value)
    except ValueError as exc:
        msg = f"{value} is not a valid number of processors"
        raise click.BadParameter(msg) from exc

    return min(number_of_processors, AVAILABLE_PROCESSORS)


def smart_labels(labels: list[str]) -> list[str]:
    inferred = [Path(label).stem for label in labels]
    if len(inferred) != len(set(inferred)):
        print(
            "Labels inferred from file names are not unique. "
            "Please be aware that in case of overlapping barcodes the counts will be "
            "merged."
        )
    return inferred


def version_callback(ctx: click.Context, value: bool) -> None:
    """Print ``<command-path> <version>`` and exit.

    Used as the callback of the shared ``--version`` option (``OTHER_OPTS["version"]``).
    The command path is reconstructed from the click context, so it works for the
    standalone console scripts (e.g. ``scClusterCells``) as well as for
    ``sincei scClusterCells --version``.
    """
    if value:
        names: list[str] = []
        node: click.Context | None = ctx
        while node is not None:
            if node.info_name:
                names.append(node.info_name)
            node = node.parent
        names.reverse()
        typer.echo(f"{' '.join(names)} {version_string()}")
        raise typer.Exit()


def help_callback(ctx: click.Context, value: bool) -> None:
    """Print the command help and exit.

    Used as the callback of the shared ``--help`` option (``OTHER_OPTS["help"]``) so
    that ``-h/--help`` shows up in the "Other options" panel instead of typer's
    default (ungrouped) help option. Commands that use it must disable typer's
    auto-generated help by passing ``context_settings={"help_option_names": []}``.
    """
    if value:
        typer.echo(ctx.get_help())
        raise typer.Exit()


def override(
    option: typer.models.OptionInfo, **changes: object
) -> typer.models.OptionInfo:
    """Return a copy of a shared ``typer.Option`` with some attributes replaced.

    Used to reuse a centrally defined option (flags, help, metavar) while giving it
    a command-specific panel or help text, e.g.
    ``override(READ_OPTS["min_mapping_quality"], rich_help_panel=_BARCODE)``.

    A command-specific *default* needs no override: in ``Annotated`` style the
    default lives in the function signature, so ``bin_size: BinSize = 100`` is all
    it takes.
    """
    new_option = copy.deepcopy(option)
    for key, value in changes.items():
        setattr(new_option, key, value)
    return new_option


# ---------------------------------------------------------------------------
# Enums for choice-style options
# ---------------------------------------------------------------------------
class PlotFileFormat(str, Enum):
    png = "png"
    jpg = "jpg"
    svg = "svg"
    pdf = "pdf"


class DuplicateFilter(str, Enum):
    start_bc = "start_bc"
    start_bc_umi = "start_bc_umi"
    start_end_bc = "start_end_bc"
    start_end_bc_umi = "start_end_bc_umi"


class FilterRNAStrand(str, Enum):
    forward = "forward"
    reverse = "reverse"


class OutFileFormat(str, Enum):
    bigwig = "bigwig"
    bedgraph = "bedgraph"


class NormalizeUsing(str, Enum):
    cpm = "CPM"
    rpkm = "RPKM"
    frequency = "Frequency"
    mean = "Mean"
    none = "None"


class Compression(str, Enum):
    none = "none"
    gzip = "gzip"


class ExportFormat(str, Enum):
    bm = "bm"
    mtx = "mtx"
    loom = "loom"


class OverlapPolicy(str, Enum):
    partial = "partial"
    all = "all"
    none = "none"


class GLMPCAFamily(str, Enum):
    poisson = "poisson"
    nb = "nb"
    mult = "mult"
    bern = "bern"


class DimRed(str, Enum):
    LSA = "LSA"
    LDA = "LDA"
    logPCA = "logPCA"
    glmPCA = "glmPCA"


class CombineMethod(str, Enum):
    # Hyphens are not valid in identifiers, so the member names differ from the
    # values the user types.
    multi_sample = "multi-sample"
    multi_modal = "multi-modal"


# ---------------------------------------------------------------------------
# Panel name constants (used by all option dicts)
# ---------------------------------------------------------------------------
_IO = "Input / Output options"
_BAM = "BAM options"
_FILTER = "Filter options"
_READ = "Read options"
_PLOT = "Plot options"
_GTF = "GTF / GFF options"
_OTHER = "Other options"


# ---------------------------------------------------------------------------
# Reusable option definitions
# ---------------------------------------------------------------------------
INPUT_OUTPUT_OPTS: dict[str, typer.models.OptionInfo] = {
    "h5ad_file": typer.Option(
        "-i",
        "--input",
        metavar=".h5ad",
        rich_help_panel=_IO,
        help="sincei-generated input file in .h5ad format.",
    ),
    "h5ad_files": typer.Option(
        "-i",
        "--input",
        metavar=".h5ad",
        rich_help_panel=_IO,
        help="List of sincei-generated input .h5ad files separated by spaces.",
    ),
    "h5mu_file": typer.Option(
        "-i",
        "--input",
        metavar=".h5mu",
        rich_help_panel=_IO,
        help="sincei-generated input file in .h5mu format.",
    ),
    "bam_file": typer.Option(
        "-b",
        "--bamfile",
        metavar=".bam",
        rich_help_panel=_IO,
        help="Indexed BAM file.",
    ),
    "bam_files": typer.Option(
        "-b",
        "--bamfiles",
        metavar=".bam",
        rich_help_panel=_IO,
        help="List of input indexed BAM files separated by spaces.",
    ),
    "barcodes": typer.Option(
        "-bc",
        "--barcodes",
        metavar=".txt",
        rich_help_panel=_IO,
        help=(
            "A single-column file containing the cell barcode whitelist, one barcode"
            "per line."
        ),
    ),
    "whitelist": typer.Option(
        "-w",
        "--whitelist",
        metavar=".txt",
        rich_help_panel=_IO,
        help="A single-column file containing the whitelist of barcodes to be used.",
    ),
    "out_file": typer.Option(
        "-o",
        "--outFile",
        metavar="PATH",
        rich_help_panel=_IO,
        help=(
            "The file to write results to. This can be a .tsv file with the results, "
            "an image for a plot, or an updated .h5ad or .h5mu file with the result of "
            "the requested operation, depending on the tool."
        ),
    ),
    "out_prefix": typer.Option(
        "-o",
        "--outFilePrefix",
        metavar="PATH",
        rich_help_panel=_IO,
        help=(
            "Prefix for output file names. One file per group is written, named "
            "`<prefix>_<group>.bw` or `<prefix>_<group>.bedgraph`."
        ),
    ),
    "group_info": typer.Option(
        "-gi",
        "--groupInfo",
        metavar=".tsv",
        rich_help_panel=_IO,
        help=(
            "A 4-column tsv file with cell grouping information in the format: "
            "`sample::barcode, UMAP1, UMAP2, group` (like the output from "
            "scClusterCells) or 3-column tsv file with format: `sample, barcode, "
            "group`. Coverages will be computed per group."
        ),
    ),
    "bed_gff_file": typer.Option(
        "-bed",
        "-gff",
        "-gtf",
        "--features",
        metavar=".bed",
        rich_help_panel=_IO,
        help=(
            "BED, GFF or GTF files to limit the coverage analysis to the regions in "
            "them."
        ),
    ),
    "region": typer.Option(
        "-r",
        "--region",
        callback=normalize_region,
        metavar="chr:start-end",
        rich_help_panel=_IO,
        help=(
            "Region of the genome to limit the operation to. Format as "
            "chr:start:end, for example ``--region chr10`` or ``--region "
            "chr10:456700:891000``.\n\n"
            "Coordinates are [bold yellow]0-based, half-open[/bold yellow], as "
            "in the BED format."
        ),
    ),
}


BAM_OPTS: dict[str, typer.models.OptionInfo] = {
    "cell_tag": typer.Option(
        "-ct",
        "--cellTag",
        metavar="XX",
        rich_help_panel=_BAM,
        help="Name of the BAM tag from which to extract barcodes.",
    ),
    "umi_tag": typer.Option(
        "-ut",
        "--umiTag",
        metavar="XX",
        rich_help_panel=_BAM,
        help=(
            "Name of the BAM tag holding the UMI, used by ``--duplicateFilter`` "
            "([bold yellow]start_bc_umi[/bold yellow] and "
            "[bold yellow]start_end_bc_umi[/bold yellow]). "
        ),
    ),
    "group_tag": typer.Option(
        "-gt",
        "--groupTag",
        metavar="XX",
        rich_help_panel=_BAM,
        help=(
            "Name of the BAM tag holding each read's sample of origin, for a BAM "
            "merged from several samples (usually ``RG`` or ``SM``).\n\n"
            "In a merged file a cell barcode may no longer be unique, since the "
            "same barcode in two source samples is two different cells. This tag, "
            "together with ``--cellTag``, is what identifies a cell.\n\n"
            "Requires a single input BAM, which must declare its samples as ``@RG`` "
            "header lines."
        ),
    ),
    "labels": typer.Option(
        "-l",
        "--labels",
        metavar="TEXT",
        rich_help_panel=_BAM,
        help=(
            "User defined labels instead of default labels from file names. Multiple "
            "labels must be separated by a space, e.g. "
            "``--labels sample1 sample2 sample3``."
        ),
    ),
    "smart_labels": typer.Option(
        "--smartLabels",
        rich_help_panel=_BAM,
        help=(
            "Instead of manually specifying labels for the input BAM files, use the "
            "file name after removing the path and extension."
        ),
    ),
    "blacklist": typer.Option(
        "-bl",
        "--blacklist",
        metavar=".bed",
        rich_help_panel=_BAM,
        help=(
            "A BED or GTF file containing regions that should be excluded from the "
            "analyses. A read is rejected if at least half of it lies inside a "
            "blacklist entry."
        ),
    ),
    "chr_to_skip": typer.Option(
        "--chrToSkip",
        metavar="CHR",
        rich_help_panel=_BAM,
        help=(
            "A space separated list of chromosomes to exclude, e.g. "
            "``--chr-to-skip chrM chrX``. Useful for skipping mitochondrial, sex "
            "chromosomes or unplaced contigs."
        ),
    ),
    "bin_size": typer.Option(
        "-bs",
        "--binSize",
        metavar="INT",
        rich_help_panel=_BAM,
        help="Size of the bins, in bases, to calculate coverage.",
    ),
    "distance_between_bins": typer.Option(
        "--distanceBetweenBins",
        metavar="INT",
        rich_help_panel=_BAM,
        help=(
            "The gap length, in bases, between bins for calculating coverage. Larger"
            "values can be used to sample the genome for input files with high "
            "coverage."
        ),
    ),
}


FILTER_OPTS: dict[str, typer.models.OptionInfo] = {
    "duplicate_filter": typer.Option(
        "--duplicateFilter",
        metavar="FILTER",
        rich_help_panel=_FILTER,
        help=(
            "How to filter for duplicates? Different combinations "
            "(using start/end/umi) "
            "are possible. Read start position and read barcode are always considered. "
            "Default (none) considers all reads as passing the filter. Note that for "
            "paired end data both reads in the fragment are considered (and kept) so, "
            "to keep only read1 combine this with ``--sam-flag-include``.\n\n"
            "One of:"
            "[bold yellow]start_bc[/bold yellow], "
            "[bold yellow]start_bc_umi[/bold yellow], "
            "[bold yellow]start_end_bc[/bold yellow], "
            "[bold yellow]start_end_bc_umi[/bold yellow]."
        ),
    ),
    "motif_filter": typer.Option(
        "-m",
        "--motifFilter",
        metavar="STR",
        rich_help_panel=_FILTER,
        help=(
            "Check whether a given motif is present in the read and the corresponding "
            "reference genome. This checks for the motif at the 5-end of the read and "
            "at the 5-overhang in the genome, which is useful in identifying reads "
            "properly cut by a restriction-enzyme or MNase. For example, to search for "
            'an "A" at the 5\'-end of the read and "TA" at the 5\'-overhang, use '
            "``-m 'A,TA'``. Reads not containing the motif are filtered out."
        ),
    ),
    "genome_2bit": typer.Option(
        "-g",
        "--genome2bit",
        metavar="file.2bit",
        rich_help_panel=_FILTER,
        help=(
            "If ``--motif-filter`` is provided, please also provide the genome "
            "sequence in 2bit format."
        ),
    ),
    "gc_content_filter": typer.Option(
        "-gc",
        "--GCcontentFilter",
        metavar="STR",
        rich_help_panel=_FILTER,
        help=(
            "Check whether the GC content of the read falls within the provided range. "
            "Input format must be '<low>,<high>' where <low> and <high> are the lower "
            "and upper bounds in the fraction of GC (e.g. '0.1,0.9'). Reads whose GC "
            "content falls outside the range are filtered out."
        ),
    ),
    "min_aligned_fraction": typer.Option(
        "--minAlignedFraction",
        metavar="FLOAT",
        rich_help_panel=_FILTER,
        help=(
            "Minimum fraction of the read that should be aligned to be counted. This "
            "includes mismatches tolerated by the aligners, but excludes "
            "InDels/clippings."
        ),
    ),
}


READ_OPTS: dict[str, typer.models.OptionInfo] = {
    "min_mapping_quality": typer.Option(
        "-mq",
        "--minMappingQuality",
        metavar="INT",
        rich_help_panel=_READ,
        help=(
            "If set, only reads that have a mapping quality score of at least this are "
            "considered. A read whose MAPQ is unavailable (255) has no score to "
            "compare and is dropped, so passing 0 keeps every scored read and drops "
            "the unscored ones rather than turning the filter off."
        ),
    ),
    "sam_flag_include": typer.Option(
        "--samFlagInclude",
        metavar="INT",
        rich_help_panel=_READ,
        help=(
            "Include reads based on SAM flag. For example, to get only reads that are "
            "the first mate, use a flag of 64. This is useful to count properly paired "
            "reads only once."
        ),
    ),
    "sam_flag_exclude": typer.Option(
        "--samFlagExclude",
        metavar="INT",
        rich_help_panel=_READ,
        help=(
            "Exclude reads based on the SAM flag. For example, to get only reads that "
            "map to the forward strand, use ``--sam-flag-exclude 16``, where 16 is the "
            "SAM flag for reads that map to the reverse strand."
        ),
    ),
    "min_fragment_length": typer.Option(
        "--minFragmentLength",
        metavar="INT",
        rich_help_panel=_READ,
        help=(
            "The minimum fragment length needed for read/pair inclusion. Useful in "
            "ATAC-seq experiments for filtering mono- or di-nucleosome fragments."
        ),
    ),
    "max_fragment_length": typer.Option(
        "--maxFragmentLength",
        metavar="INT",
        rich_help_panel=_READ,
        help="The maximum fragment length accepted for read/pair inclusion.",
    ),
    "filter_rna_strand": typer.Option(
        "--filterRNAstrand",
        metavar="STRAND",
        rich_help_panel=_READ,
        help=(
            "Selects RNA-seq reads (single-end or paired-end) originating from genes "
            "on the given strand. This assumes a standard dUTP-based library "
            "preparation (that is, ``--filter-rna-strand forward`` keeps minus-strand "
            "reads, which originally came from genes on the forward strand using a "
            "dUTP-based method). Consider using ``--sam-flag-exclude`` instead for "
            "filtering by strand in other contexts.\n\n"
            "One of: "
            "[bold yellow]forward[/bold yellow], "
            "[bold yellow]reverse[/bold yellow]."
        ),
    ),
    "extend_reads": typer.Option(
        "-e",
        "--extendReads",
        metavar="INT",
        rich_help_panel=_READ,
        help=(
            "Extend reads to a fragment size (in bases). Pass a value to use it "
            "directly, or pass the flag with no value (``--extendReads``) to estimate "
            "the fragment size from the data (the median template length of properly "
            "paired reads; errors on single-end data). If omitted, reads are not "
            "extended.\n\n"
            "*NOTE*: generally NOT recommended for spliced-read data such as RNA-seq, "
            "as it would extend reads over skipped regions.\n\n"
            "*Single-end*: the value is the final fragment length; reads that already "
            "exceed it are left unchanged.\n\n"
            "*Paired-end*: reads with mates are extended to match the fragment size "
            "defined by the two read mates (mates that map too far apart, >4x the "
            "fragment size, or to different chromosomes are treated as single-end)."
        ),
    ),
    "center_reads": typer.Option(
        "--centerReads",
        rich_help_panel=_READ,
        help=(
            "Center reads with respect to the fragment length. For paired-end data the "
            "read is centered at the fragment length defined by the two ends of the "
            "fragment; for single-end data the given fragment length is used. Useful "
            "to get a sharper signal around enriched regions."
            "*NOTE*: generally NOT recommended for spliced-read data such as RNA-seq, "
            "as it would make reads cover skipped regions.\n\n"
        ),
    ),
}


OTHER_OPTS: dict[str, typer.models.OptionInfo] = {
    "number_of_processors": typer.Option(
        "-p",
        "--numberOfProcessors",
        parser=normalize_processors,
        metavar="INT",
        show_default="max",
        rich_help_panel=_OTHER,
        help=(
            'Number of processors to use. You can also type "max/2" to use half the '
            'maximum number of processors or "max" to use all available processors. '
            '(Default: "max")'
        ),
    ),
    "verbose": typer.Option(
        "-v",
        "--verbose",
        rich_help_panel=_OTHER,
        help="Set to see processing messages.",
    ),
    "version": typer.Option(
        "-V",
        "--version",
        is_eager=True,
        expose_value=False,
        callback=version_callback,
        rich_help_panel=_OTHER,
        help="Print the program version and exit.",
    ),
    "help": typer.Option(
        "-h",
        "--help",
        is_eager=True,
        expose_value=False,
        callback=help_callback,
        rich_help_panel=_OTHER,
        help="Show this message and exit.",
    ),
}


PLOT_OPTS: dict[str, typer.models.OptionInfo] = {
    "plot_width": typer.Option(
        "--plotWidth",
        metavar="FLOAT",
        rich_help_panel=_PLOT,
        help="Output plot width (in cm).",
    ),
    "plot_height": typer.Option(
        "--plotHeight",
        metavar="FLOAT",
        rich_help_panel=_PLOT,
        help="Output plot height (in cm).",
    ),
    "plot_file_format": typer.Option(
        "--plotFileFormat",
        metavar="FORMAT",
        rich_help_panel=_PLOT,
        help=(
            "Image format type. If given, this option will be used to save the plot "
            "file.\n\n"
            "One of: [bold yellow]png[/bold yellow], [bold yellow]jpg[/bold yellow], "
            "[bold yellow]svg[/bold yellow], [bold yellow]pdf[/bold yellow]."
        ),
    ),
}


GTF_GFF_OPTS: dict[str, typer.models.OptionInfo] = {
    "metagene": typer.Option(
        "-m",
        "--metagene",
        rich_help_panel=_GTF,
        help=(
            "When a GFF/GTF file is used to provide regions, count reads only on the "
            "combined exons of a transcript or gene rather than on the genomic "
            "interval defined by the 5-prime and 3-prime transcript bound. Exons  are"
            "the records whose column-3 type matches ``--exonID``; what they are "
            "grouped into is set by ``--transcriptIDtag``: one feature per transcript "
            "(the default) or one per gene (``--transcriptIDtag gene_id``). "
            "Ignored for BED inputs."
        ),
    ),
    "transcript_id": typer.Option(
        "--transcriptID",
        metavar="STR",
        rich_help_panel=_GTF,
        help=(
            "The column-3 feature type(s) processed as a region (transcript). May be "
            "given more than once, e.g. ``--transcriptID mRNA --transcriptID lnc_RNA``."
            " (Default: every transcript type in the file, regardless of "
            "biotype. GTF and GENCODE GFF3 type them all ``transcript``; an "
            "Ensembl-style GFF3 names them by biotype instead: mRNA, lnc_RNA, snoRNA, "
            "etc., which are read out of the file itself.)"
        ),
    ),
    "exon_id": typer.Option(
        "--exonID",
        metavar="STR",
        rich_help_panel=_GTF,
        help=(
            "When a GTF/GFF file is used to provide regions, entries with this value "
            "as their feature (column 3) are treated as exons. May be given more than "
            'once, e.g. ``--exonID exon --exonID CDS``. "CDS" is another common value. '
            "NOTE: only used in metagene mode. (Default: ``exon``, which all three "
            "flavours agree on.)"
        ),
    ),
    "transcript_id_tag": typer.Option(
        "--transcriptIDtag",
        metavar="STR",
        rich_help_panel=_GTF,
        help=(
            "The column-9 attribute key whose value names each region. For GTF this is "
            "stored as a key:value pair in the last column. Set this to use a "
            "different identifier. (Default: ``transcript_id`` for GTF, ``ID`` for "
            "GFF3.)\n\n"
            "In ``--metagene`` mode this key instead selects what exons are grouped "
            "by, and so what one feature is. The default groups per transcript "
            "(``transcript_id`` for GTF, ``Parent`` for GFF3, both of which every "
            "flavour carries). Pass ``--transcriptIDtag gene_id`` to group per gene "
            "instead. Note that an Ensembl-style GFF3 does not carry a gene id on "
            "its exons, so grouping by gene gives an error."
        ),
    ),
}
