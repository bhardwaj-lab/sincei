"""Shared CLI building blocks for the sincei argparse apps.

This module is the argparse counterpart of ``cli_typer/_common_args.py`` and is
kept flag-for-flag identical to it: same option names, same defaults, same
choices and the same help text, grouped into the same panels.

Each factory returns a parent ``ArgumentParser`` (``add_help=False``) that the
command modules combine via ``parents=[...]``.

Note on selecting options: every factory takes an explicit ``include`` set of
option names.  An option that is not included is genuinely absent from the
parser.  (The previous implementation passed ``argparse.SUPPRESS`` as the *help
text*, which only hid the option from ``--help`` while still accepting it on the
command line.)
"""

from __future__ import annotations

import argparse
import os
from typing import TYPE_CHECKING

from sincei import _sincei as internal

from ._parsers import DUPLICATE_FILTER_CHOICES

if TYPE_CHECKING:
    from collections.abc import Iterable

# ---------------------------------------------------------------------------
# Choice-style options (the argparse counterpart of the typer enums)
# ---------------------------------------------------------------------------
PLOT_FILE_FORMATS = ["png", "jpg", "svg", "pdf"]
FILTER_RNA_STRANDS = ["forward", "reverse"]
OUT_FILE_FORMATS = ["bigwig", "bedgraph"]
NORMALIZE_USING = ["CPM", "RPKM", "Frequency", "Mean", "None"]
COMPRESSION_METHODS = ["none", "gzip"]
EXPORT_FORMATS = ["bm", "mtx", "loom"]
OVERLAP_POLICIES = ["partial", "all", "none"]
GLMPCA_FAMILIES = ["poisson", "nb", "mult", "bern"]
DIM_RED_METHODS = ["LSA", "LDA", "logPCA", "glmPCA"]

# Panel (argument-group) titles, matching the typer help panels.
_IO = "Input / Output options"
_BAM = "BAM options"
_FILTER = "Filter options"
_READ = "Read options"
_PLOT = "Plot options"
_GTF = "GTF / GFF options"
_OTHER = "Other options"


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def available_processors() -> int:
    """Usable processor count (``os.sched_getaffinity`` where it exists)."""
    try:
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return os.cpu_count() or 1


AVAILABLE_PROCESSORS = available_processors()


def genomic_region(string: str) -> str | None:
    """``--region`` type: normalize ``chr:start-end`` into ``chr:start:end``."""
    region = "".join(string.split())
    if region == "":
        return None
    region = region.translate({ord(c): None for c in ",;|!{}()"}).replace("-", ":")
    if len(region) == 0:
        msg = f"{string} is not a valid region"
        raise argparse.ArgumentTypeError(msg)
    return region


def number_of_processors(string: str) -> int:
    """``--numberOfProcessors`` type: accept ``max``, ``max/2`` or an integer."""
    if string == "max/2":
        return int(AVAILABLE_PROCESSORS * 0.5)
    if string == "max":
        return AVAILABLE_PROCESSORS
    try:
        value = int(string)
    except ValueError:
        msg = f"{string} is not a valid number of processors"
        raise argparse.ArgumentTypeError(msg) from None
    return min(value, AVAILABLE_PROCESSORS)


def _wanted(name: str, include: Iterable[str] | None) -> bool:
    return include is None or name in include


def _default(name: str, defaults: dict[str, object] | None, fallback: object) -> object:
    if defaults and name in defaults:
        return defaults[name]
    return fallback


# ---------------------------------------------------------------------------
# Input / Output options
# ---------------------------------------------------------------------------
def input_output_options(
    include: Iterable[str], defaults: dict[str, object] | None = None
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(_IO)

    if _wanted("h5ad_file", include):
        group.add_argument(
            "-i",
            "--input",
            metavar=".h5ad",
            required=True,
            help="sincei-generated input file in .h5ad format.",
        )
    if _wanted("h5ad_files", include):
        group.add_argument(
            "-i",
            "--input",
            metavar=".h5ad",
            action="extend",
            nargs="+",
            required=True,
            help="List of sincei-generated input .h5ad files separated by spaces.",
        )
    if _wanted("h5mu_file", include):
        group.add_argument(
            "-i",
            "--input",
            metavar=".h5mu",
            required=True,
            help="sincei-generated input file in .h5mu format.",
        )
    if _wanted("bam_file", include):
        group.add_argument(
            "-b",
            "--bamfile",
            metavar=".bam",
            required=True,
            help="Indexed BAM file.",
        )
    if _wanted("bam_files", include):
        group.add_argument(
            "-b",
            "--bamfiles",
            metavar=".bam",
            action="extend",
            nargs="+",
            required=True,
            help="List of input indexed BAM files separated by spaces.",
        )
    if _wanted("barcodes", include):
        group.add_argument(
            "-bc",
            "--barcodes",
            metavar=".txt",
            required=True,
            help="A single-column file containing the cell barcode whitelist, one "
            "barcode per line.",
        )
    if _wanted("whitelist", include):
        group.add_argument(
            "-w",
            "--whitelist",
            metavar=".txt",
            default=None,
            help="A single-column file containing the whitelist of barcodes to be "
            "used.",
        )
    if _wanted("out_file", include):
        group.add_argument(
            "-o",
            "--outFile",
            dest="outFile",
            metavar=_default("out_file_metavar", defaults, "PATH"),
            required=True,
            help=_default(
                "out_file_help",
                defaults,
                "The file to write results to. This can be a .tsv file with the "
                "results, an image for a plot, or an updated .h5ad or .h5mu file with "
                "the result of the requested operation, depending on the tool.",
            ),
        )
    if _wanted("out_prefix", include):
        group.add_argument(
            "-o",
            "--outFilePrefix",
            dest="outFilePrefix",
            metavar="PATH",
            required=True,
            help="Prefix for output file names.",
        )
    if _wanted("group_info", include):
        group.add_argument(
            "-gi",
            "--groupInfo",
            dest="groupInfo",
            metavar=".tsv",
            required=True,
            help="A 4-column tsv file with cell grouping information in the format: "
            "`sample::barcode, UMAP1, UMAP2, group` (like the output from "
            "scClusterCells) or 3-column tsv file with format: `sample, barcode, "
            "group`. Coverages will be computed per group.",
        )
    if _wanted("region", include):
        group.add_argument(
            "-r",
            "--region",
            metavar="chr:start-end",
            type=genomic_region,
            required=bool(_default("region_required", defaults, False)),
            default=None,
            help="Region of the genome to limit the operation to - this is useful "
            "when testing parameters to reduce the computing time. The format is "
            "chr:start:end, for example ``--region chr10`` or "
            "``--region chr10:456700:891000``.",
        )
    if _wanted("compression", include):
        group.add_argument(
            "--compression",
            metavar="METHOD",
            choices=COMPRESSION_METHODS,
            default="none",
            help="HDF5 compression for the output .h5ad datasets.\n\n"
            "none: no compression (default; fastest, larger files).\n\n"
            "gzip: deflate compression (smaller filesize).",
        )
        group.add_argument(
            "--compressionLevel",
            dest="compressionLevel",
            metavar="INT",
            type=int,
            default=4,
            help="Compression level (0-9) when --compression is gzip. Ignored for "
            "'none'.",
        )
    if _wanted("bed", include):
        group.add_argument(
            "--bed",
            metavar=".bed/.gtf/.gff",
            required=True,
            help="BED/GTF/GFF files to limit the coverage analysis to the regions in "
            "them.",
        )
    return parser


# ---------------------------------------------------------------------------
# BAM options
# ---------------------------------------------------------------------------
def bam_options(
    include: Iterable[str],
    defaults: dict[str, object] | None = None,
    panel: str = _BAM,
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(panel)

    if _wanted("cell_tag", include):
        group.add_argument(
            "-ct",
            "--cellTag",
            dest="cellTag",
            metavar="XX",
            default="BC",
            help="Name of the BAM tag from which to extract barcodes.",
        )
    if _wanted("group_tag", include):
        group.add_argument(
            "-gt",
            "--groupTag",
            dest="groupTag",
            metavar="XX",
            default=None,
            help="In case of a grouped BAM file, such as the one containing Read Group "
            "(``RG``) or Sample (``SM``) tag, it is possible to group the reads using "
            "the provided ``--grouptag`` argument. NOTE: in case of such input, please "
            "ensure that the ``--labels`` argument indicates the expected group labels "
            "contained in the BAM files. The ``--grouptag`` along with the "
            "``--celltag`` is then used to identify unique samples (cells) from the "
            "input.",
        )
    if _wanted("labels", include):
        group.add_argument(
            "-l",
            "--labels",
            metavar="TEXT",
            action="extend",
            nargs="+",
            required=bool(_default("labels_required", defaults, False)),
            default=None,
            help="User defined labels instead of default labels from file names. "
            "Multiple labels must be separated by a space, e.g. "
            "``--labels sample1 sample2 sample3``.",
        )
    if _wanted("smart_labels", include):
        group.add_argument(
            "--smartLabels",
            dest="smartLabels",
            action="store_true",
            help="Instead of manually specifying labels for the input BAM files, use "
            "the file name after removing the path and extension.",
        )
    if _wanted("blacklist", include):
        group.add_argument(
            "-bl",
            "--blacklist",
            metavar=".bed",
            action="extend",
            nargs="+",
            default=None,
            help="A BED or GTF file containing regions that should be excluded from "
            "all analyses. Currently this works by rejecting genomic chunks that "
            "happen to overlap an entry, so a read partially overlapping a blacklisted "
            "region (or a fragment spanning it) might still be considered. Adjust the "
            "effective genome size if relevant.",
        )
    if _wanted("chr_to_skip", include):
        group.add_argument(
            "--chrToSkip",
            dest="chrToSkip",
            metavar="CHR",
            action="extend",
            nargs="+",
            default=None,
            help="A space separated list of chromosomes to exclude, e.g. "
            "``--chrToSkip chrM chrX``. Useful for skipping mitochondrial, sex "
            "chromosomes or unplaced contigs.",
        )
    if _wanted("bin_size", include):
        group.add_argument(
            "-bs",
            "--binSize",
            dest="binSize",
            metavar="INT",
            type=int,
            default=_default("bin_size", defaults, 10000),
            help="Size of the bins, in bases, to calculate coverage.",
        )
    if _wanted("distance_between_bins", include):
        group.add_argument(
            "--distanceBetweenBins",
            dest="distanceBetweenBins",
            metavar="INT",
            type=int,
            default=_default("distance_between_bins", defaults, 0),
            help="The gap length, in bases, between bins for calculating coverage. "
            "Larger values can be used to sample the genome for input files with high "
            "coverage.",
        )
    return parser


# ---------------------------------------------------------------------------
# Filter options
# ---------------------------------------------------------------------------
def filter_options(include: Iterable[str]) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(_FILTER)

    if _wanted("duplicate_filter", include):
        group.add_argument(
            "--duplicateFilter",
            dest="duplicateFilter",
            metavar="FILTER",
            choices=DUPLICATE_FILTER_CHOICES,
            default=None,
            help="How to filter for duplicates? Different combinations (using "
            "start/end/umi) are possible. Read start position and read barcode are "
            "always considered. Default (none) considers all reads as passing the "
            "filter. Note that for paired end data both reads in the fragment are "
            "considered (and kept) so, to keep only read1 combine this with "
            "``--samFlagInclude``.",
        )
    if _wanted("motif_filter", include):
        group.add_argument(
            "-m",
            "--motifFilter",
            dest="motifFilter",
            metavar="STR",
            action="extend",
            nargs="+",
            default=None,
            help="Check whether a given motif is present in the read and the "
            "corresponding reference genome. This checks for the motif at the 5-end of "
            "the read and at the 5-overhang in the genome, which is useful in "
            "identifying reads properly cut by a restriction-enzyme or MNase. For "
            'example, to search for an "A" at the 5\'-end of the read and "TA" at '
            "the 5'-overhang, use ``-m 'A,TA'``. Reads not containing the motif are "
            "filtered out.",
        )
    if _wanted("genome_2bit", include):
        group.add_argument(
            "-g",
            "--genome2bit",
            dest="genome2bit",
            metavar="file.2bit",
            default=None,
            help="If ``--motifFilter`` is provided, please also provide the genome "
            "sequence in 2bit format.",
        )
    if _wanted("gc_content_filter", include):
        group.add_argument(
            "-gc",
            "--GCcontentFilter",
            dest="GCcontentFilter",
            metavar="STR",
            default=None,
            help="Check whether the GC content of the read falls within the provided "
            "range. Input format must be '<low>,<high>' where <low> and <high> are "
            "the lower and upper bounds in the fraction of GC (e.g. '0.1,0.9'). Reads "
            "whose GC content falls outside the range are filtered out.",
        )
    if _wanted("min_aligned_fraction", include):
        group.add_argument(
            "--minAlignedFraction",
            dest="minAlignedFraction",
            metavar="FLOAT",
            type=float,
            default=None,
            help="Minimum fraction of the read that should be aligned to be counted. "
            "This includes mismatches tolerated by the aligners, but excludes "
            "InDels/clippings.",
        )
    return parser


# ---------------------------------------------------------------------------
# Read options
# ---------------------------------------------------------------------------
def read_options(include: Iterable[str]) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(_READ)

    if _wanted("min_mapping_quality", include):
        group.add_argument(
            "-mq",
            "--minMappingQuality",
            dest="minMappingQuality",
            metavar="INT",
            type=int,
            default=None,
            help="If set, only reads that have a mapping quality score of at least "
            "this are considered.",
        )
    if _wanted("sam_flag_include", include):
        group.add_argument(
            "--samFlagInclude",
            dest="samFlagInclude",
            metavar="INT",
            type=int,
            default=None,
            help="Include reads based on SAM flag. For example, to get only reads that "
            "are the first mate, use a flag of 64. This is useful to count properly "
            "paired reads only once.",
        )
    if _wanted("sam_flag_exclude", include):
        group.add_argument(
            "--samFlagExclude",
            dest="samFlagExclude",
            metavar="INT",
            type=int,
            default=None,
            help="Exclude reads based on the SAM flag. For example, to get only reads "
            "that map to the forward strand, use ``--samFlagExclude 16``, where 16 is "
            "the SAM flag for reads that map to the reverse strand.",
        )
    if _wanted("min_fragment_length", include):
        group.add_argument(
            "--minFragmentLength",
            dest="minFragmentLength",
            metavar="INT",
            type=int,
            default=0,
            help="The minimum fragment length needed for read/pair inclusion. Useful "
            "in ATAC-seq experiments for filtering mono- or di-nucleosome fragments.",
        )
    if _wanted("max_fragment_length", include):
        group.add_argument(
            "--maxFragmentLength",
            dest="maxFragmentLength",
            metavar="INT",
            type=int,
            default=0,
            help="The maximum fragment length accepted for read/pair inclusion.",
        )
    if _wanted("filter_rna_strand", include):
        group.add_argument(
            "--filterRNAstrand",
            dest="filterRNAstrand",
            metavar="STRAND",
            choices=FILTER_RNA_STRANDS,
            default=None,
            help="Selects RNA-seq reads (single-end or paired-end) originating from "
            "genes on the given strand. This assumes a standard dUTP-based library "
            "preparation (that is, ``--filterRNAstrand forward`` keeps minus-strand "
            "reads, which originally came from genes on the forward strand using a "
            "dUTP-based method). Consider using ``--samFlagExclude`` instead for "
            "filtering by strand in other contexts.",
        )
    if _wanted("extend_reads", include):
        group.add_argument(
            "-e",
            "--extendReads",
            dest="extendReads",
            metavar="INT",
            type=int,
            nargs="?",
            const=0,
            default=None,
            help="Extend reads to a fragment size (in bases). Pass a value to use it "
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
            "fragment size, or to different chromosomes are treated as single-end).",
        )
    if _wanted("center_reads", include):
        group.add_argument(
            "--centerReads",
            dest="centerReads",
            action="store_true",
            help="Center reads with respect to the fragment length. For paired-end "
            "data the read is centered at the fragment length defined by the two ends "
            "of the fragment; for single-end data the given fragment length is used. "
            "Useful to get a sharper signal around enriched regions.",
        )
    return parser


# ---------------------------------------------------------------------------
# GTF / GFF options
# ---------------------------------------------------------------------------
def gtf_gff_options(
    include: Iterable[str], metagene_flags: tuple[str, ...] = ("--metagene",)
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(_GTF)

    if _wanted("transcript_id", include):
        group.add_argument(
            "--transcriptID",
            dest="transcriptID",
            metavar="STR",
            action="extend",
            nargs="+",
            default=None,
            help="The column-3 feature type(s) processed as a region (transcript). May "
            "be given more than once, e.g. ``--transcriptID mRNA --transcriptID "
            "lnc_RNA``. (Default: every transcript type in the file, regardless of "
            "biotype. GTF and GENCODE GFF3 type them all ``transcript``; an "
            "Ensembl-style GFF3 names them by biotype instead: mRNA, lnc_RNA, snoRNA, "
            "etc., which are read out of the file itself.)",
        )
    if _wanted("exon_id", include):
        group.add_argument(
            "--exonID",
            dest="exonID",
            metavar="STR",
            action="extend",
            nargs="+",
            default=None,
            help="When a GTF/GFF file is used to provide regions, entries with this "
            "value as their feature (column 3) are treated as exons. May be given more "
            'than once, e.g. ``--exonID exon --exonID CDS``. "CDS" is another common '
            "value. NOTE: only used in metagene mode. (Default: ``exon``, which all "
            "three flavours agree on.)",
        )
    if _wanted("transcript_id_tag", include):
        group.add_argument(
            "--transcriptIDtag",
            dest="transcriptIDtag",
            metavar="STR",
            default=None,
            help="The column-9 attribute key whose value names each region. For GTF "
            "this is stored as a key:value pair in the last column. Set this to use a "
            "different identifier. (Default: ``transcript_id`` for GTF, ``ID`` for "
            "GFF3.)\n\n"
            "In ``--metagene`` mode this key instead selects what exons are grouped "
            "by, and so what one feature is. The default groups per transcript "
            "(``transcript_id`` for GTF, ``Parent`` for GFF3, both of which every "
            "flavour carries). Pass ``--transcriptIDtag gene_id`` to group per gene "
            "instead -- note that an Ensembl-style GFF3 does not carry a gene id on "
            "its exons, so grouping one by gene is an error.",
        )
    if _wanted("metagene", include):
        group.add_argument(
            *metagene_flags,
            dest="metagene",
            action="store_true",
            help="When a GFF/GTF file is used to provide regions, count reads on the "
            "combined exons of a transcript or gene rather than on the genomic "
            "interval defined by the 5-prime and 3-prime most transcript bound. Exons "
            "are the records whose column-3 type matches ``--exonID``; what they are "
            "grouped into -- one feature per transcript (the default) or per gene "
            "(``--transcriptIDtag gene_id``) -- is set by ``--transcriptIDtag``. "
            "Ignored for BED inputs.",
        )
    return parser


# ---------------------------------------------------------------------------
# Plot options
# ---------------------------------------------------------------------------
def plot_options() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(_PLOT)

    group.add_argument(
        "--plotWidth",
        dest="plotWidth",
        metavar="FLOAT",
        type=float,
        default=10.0,
        help="Output plot width (in cm).",
    )
    group.add_argument(
        "--plotHeight",
        dest="plotHeight",
        metavar="FLOAT",
        type=float,
        default=10.0,
        help="Output plot height (in cm).",
    )
    group.add_argument(
        "--plotFileFormat",
        dest="plotFileFormat",
        metavar="FORMAT",
        choices=PLOT_FILE_FORMATS,
        default="png",
        help="Image format type. If given, this option will be used to save the plot "
        "file.",
    )
    return parser


# ---------------------------------------------------------------------------
# Other options
# ---------------------------------------------------------------------------
def other_options(version: bool = False) -> argparse.ArgumentParser:
    """Shared trailing options.

    ``version`` adds ``-V/--version``; only the ``sincei`` root command carries
    it, matching the typer CLI.
    """
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(_OTHER)

    group.add_argument(
        "-p",
        "--numberOfProcessors",
        dest="numberOfProcessors",
        metavar="INT",
        type=number_of_processors,
        default=AVAILABLE_PROCESSORS,
        help='Number of processors to use. You can also type "max/2" to use half the '
        'maximum number of processors or "max" to use all available processors. '
        '(Default: "max")',
    )
    group.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Set to see processing messages.",
    )
    if version:
        group.add_argument(
            "-V",
            "--version",
            action="version",
            version=f"%(prog)s {internal.version()}",
            help="Print the program version and exit.",
        )
    group.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this message and exit.",
    )
    return parser


def build_parser(
    description: str, usage: str | None, parents: list[argparse.ArgumentParser]
) -> argparse.ArgumentParser:
    """Assemble a command parser from its parent option groups."""
    return argparse.ArgumentParser(
        parents=parents,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=description,
        usage=usage,
        add_help=False,
        conflict_handler="resolve",
    )
