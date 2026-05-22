from __future__ import annotations

import os
import sys
import logging
import typer

from enum import Enum

from sincei import _sincei as internal

INPUT_OUTPUT_OPTS = {
    "h5ad_file": typer.Option(
        "-i",
        "--input",
        default=...,
        metavar=".h5ad",
        help="sincei-generated input file in .h5ad format.",
        rich_help_panel="Input/Output options",
    ),
    "h5ad_files": typer.Option(
        "-i",
        "--input",
        default=...,
        metavar=".h5ad",
        help="List of sincei-generated input .h5ad files separated by spaces.",
        rich_help_panel="Input/Output options",
    ),
    "h5mu_file": typer.Option(
        "-i",
        "--input",
        default=...,
        metavar=".h5mu",
        help="sincei-generated input file in .h5mu format.",
        rich_help_panel="Input/Output options",
    ),
    "bam_file": typer.Option(
        "-b",
        "--bam",
        default=...,
        metavar=".bam",
        help="Indexed BAM file.",
        rich_help_panel="Input/Output options",
    ),
    "bam_files": typer.Option(
        "-b",
        "--bam",
        default=...,
        metavar=".bam",
        help="List of input indexed BAM files separated by spaces.",
        rich_help_panel="Input/Output options",
    ),
    "out_file": typer.Option(
        "-o",
        "--out-file",
        default=...,
        metavar="path",
        help="The file to write results to. This can be a .tsv file with the results, an image for "
        "a plot, or an updated .h5ad or .h5mu file with the result of the requested operation, "
        "depending on the tool.",
        rich_help_panel="Input/Output options",
    ),
    "out_prefix": typer.Option(
        "-o",
        "--out-prefix",
        default=...,
        metavar="path",
        help="Prefix for output file names.",
        rich_help_panel="Input/Output options",
    ),
    "whitelist": typer.Option(
        "-wl",
        "--whitelist",
        default=None,
        metavar=".txt",
        help="A single-column file containing the whitelist of barcodes to be used.",
        rich_help_panel="Input/Output options",
    ),
    "barcodes": typer.Option(
        "-bc",
        "--barcodes",
        default=None,
        metavar=".txt",
        help="A single-column file containing the list of barcodes to be used.",
        rich_help_panel="Input/Output options",
    ),
    "bed_gff_file": typer.Option(
        "-bed",
        "-gff",
        "-gtf",
        "--intervals",
        default=None,
        metavar=".bed",
        help="BED, GFF or GTF files to limit the coverage analysis to the regions in them.",
        rich_help_panel="Input/Output options",
    ),
    "group_info": typer.Option(
        "-g",
        "--group-info",
        default=None,
        metavar=".tsv",
        help="A 4-column tsv file with cell grouping information in the format: `sample::barcode, "
        "UMAP1, UMAP2, group` (like the output from scClusterCells) or 3-column tsv file with "
        "format: `sample, barcode, group`. Coverages will be computed per group.",
        rich_help_panel="Input/Output options",
    ),
    "labels": typer.Option(
        "-l",
        "--labels",
        default=None,
        metavar="sample",
        help="User defined labels instead of default labels from file names. Multiple labels must be "
        "separated by a space, e.g. ``--labels sample1 sample2 sample3``.",
        rich_help_panel="Input/Output options",
    ),
    "smart_labels": typer.Option(
        "--smart-labels",
        default=False,
        help="Instead of manually specifying labels for the input BAM files, use the file name after "
        "removing the path and extension.",
        rich_help_panel="Input/Output options",
    ),
    "region": typer.Option(
        "-r",
        "--region",
        default=None,
        parser=parse_region,
        metavar="chr:start-end",
        help="Genomic region of the genome to limit the operation to - this is useful when testing "
        "parameters to reduce the computing time. Must be given in format chr:start:end, for example "
        "``--region chr10`` or ``--region chr10:456700:891000``.",
        rich_help_panel="Input/Output options",
    ),
}


OTHER_OPTS = {
    "number_of_processors": typer.Option(
        "-p",
        "--n-processors",
        default="max",
        parser=parse_processors,
        help='Number of processors to use. Type "max/2" to use half the maximum number of processors "'
        'or "max" to use all available processors. (Default: "max")',
        rich_help_panel="Other options",
    ),
    "help": typer.Option(
        "-h",
        "--help",
        default=None,
        help="Show this message and exit.",
        rich_help_panel="Other options",
        callback=help_option,
    ),
    "verbose": typer.Option(
        "-v",
        "--verbose",
        default=False,
        help="Set to see processing messages.",
        rich_help_panel="Other options",
    ),
    "version": typer.Option(
        "-V",
        "--version",
        default=None,
        is_eager=True,
        expose_value=False,
        help="Print the program version and exit.",
        rich_help_panel="Other options",
        callback=version_option,
    ),
}


GTF_GFF_OPTS = {
    "metagene": typer.Option(
        "-m",
        "--metagene",
        default=False,
        metavar="str",
        help="When either a BED12 or GTF file are used to provide regions, "
        "perform the computation on the merged exons, rather than using the "
        "genomic interval defined by the 5-prime and 3-prime most transcript "
        "bound (i.e., columns 2 and 3 of a BED file). If a BED3 or BED6 file is "
        "used as input, then columns 2 and 3 are used as an exon. "
        "(Default: %(default)s)",
        rich_help_panel="GTF/GFF options",
    ),
    "transcript_id": typer.Option(
        "--transcript-id",
        default="transcript",
        metavar="str",
        help="When a GTF file is used to provide regions, only  entries with this value as their "
        "feature (column 3) will be processed as transcripts. (Default: %(default)s)",
        rich_help_panel="GTF/GFF options",
    ),
    "exon_id": typer.Option(
        "--exon-id",
        default="exon",
        metavar="str",
        help="When a GTF file is used to provide regions, only entries with "
        "this value as their feature (column 3) will be processed as exons. "
        "CDS would be another common value for this. (Default: %(default)s)",
        rich_help_panel="GTF/GFF options",
    ),
    "transcript_id_tag": typer.Option(
        "--transcript-id-tag",
        default="transcript_id",
        metavar="str",
        help="Each region has an ID (e.g., ACTB) assigned to it, which for BED files is either "
        "column 4 (if it exists) or the interval bounds. For GTF files this is instead stored in "
        "the last column as a key:value pair (e.g., as 'transcript_id \"ACTB\"', for a key of "
        "transcript_id and a value of ACTB). In some cases it can be convenient to use a different "
        "identifier. To do so, set this to the desired key. (Default: %(default)s)",
        rich_help_panel="GTF/GFF options",
    ),
}

PLOT_OPTS = {
    "plot_width": typer.Option(
        "--plot-width",
        default=10.0,
        metavar="float",
        help="Output plot width (in cm). (Default: %(default)s)",
        rich_help_panel="Plot options",
    ),
    "plot_height": typer.Option(
        "--plot-height",
        default=10.0,
        metavar="float",
        help="Output plot height (in cm). (Default: %(default)s)",
        rich_help_panel="Plot options",
    ),
    "plot_file_format": typer.Option(
        "--plot-file-format",
        default=PlotFileFormat.png,
        help="Image format type. If given, this option will be used to save the plot file. (Default: %(default)s)",
        rich_help_panel="Plot options",
    ),
}

BAM_OPTS = {
    "cell_tag": typer.Option(
        "-ct",
        "--cell-tag",
        default="BC",
        metavar="XX",
        help="Name of the BAM tag from which to extract barcodes. (Default: %(default)s)",
        rich_help_panel="BAM options",
    ),
    "group_tag": typer.Option(
        "-gt",
        "--group-tag",
        default=None,
        metavar="XX",
        help="In case of a groupped BAM file, such as the one containing Read Group (``RG``) or "
        "Sample (``SM``) tag, it is possible to process group the reads using the provided "
        "``--groupTag`` argument. NOTE: In case of such input, please ensure that the ``--labels`` "
        "argument indicates the expected group labels contained in the BAM files. The "
        "``--groupTag`` along with the ``--cellTag`` is then used to identify unique samples "
        "(cells) from the input",
        rich_help_panel="BAM options",
    ),
    "labels": typer.Option(
        "-l",
        "--labels",
        default=None,
        metavar="sample",
        help="User defined labels instead of default labels from file names. Multiple labels have "
        "to be separated by a space, e.g. ``--labels sample1 sample2 sample3``.",
        rich_help_panel="BAM options",
    ),
    "smart_labels": typer.Option(
        "--smart-labels",
        default=False,
        help="Instead of manually specifying labels for the input BAM files, use the file name after "
        "removing the path and extension.",
        rich_help_panel="BAM options",
    ),
    "region": typer.Option(
        "-r",
        "--region",
        default=None,
        parser=parse_region,
        metavar="chr:start-end",
        help="Region of the genome to limit the operation to - this is useful when testing "
        "parameters to reduce the computing time. The format is chr:start:end, for example "
        "``--region chr10`` or ``--region chr10:456700:891000``.",
        rich_help_panel="BAM options",
    ),
    "blacklist": typer.Option(
        "-bl",
        "--blacklist",
        default=None,
        metavar=".bed",
        help="BED/GTF files to limit analysis to the provided regions.",
        rich_help_panel="BAM options",
    ),
    "chr_to_skip": typer.Option(
        "--chr-to-skip",
        default=None,
        metavar="chr",
        help="A space separated list of chromosomes to exclude. eg. chrM chrUn",
        rich_help_panel="BAM options",
    ),
    "binsize": typer.Option(
        "-bs",
        "--binsize",
        default=None,
        metavar="int",
        help="Bin size to use for counting reads. (Default: %(default)s)",
        rich_help_panel="BAM options",
    ),
    "distance_between_bins": typer.Option(
        "--distance-between-bins",
        default=0,
        metavar="int",
        help="Distance between bins to use for counting reads. (Default: %(default)s)",
        rich_help_panel="BAM options",
    ),
}


FILTER_OPTS = {
    "duplicate_filter": typer.Option(
        "--duplicate-filter",
        default=DuplicateFilter.none,
        help="How to filter for duplicates? Different combinations (using start/end/umi) are "
        "possible. Read start position and read barcode are always considered. Default (None) "
        "considers all reads as passing the filter.\n"
        "Note that in case of paired end data, both reads in the fragment are considered (and kept). "
        "So if you wish to keep only read1, combine this option with `--samFlagInclude`. ",
    ),
    "motif_filter": typer.Option(
        "-m",
        "--motif-filter",
        default=None,
        metavar="str",
        help="Check whether a given motif is present in the read and the corresponding reference "
        "genome. This option checks for the motif at the 5-end of the read and at the 5-overhang in "
        "the genome, which is useful in identifying reads properly cut by a restriction-enzyme or "
        'MNAse. For example, if you want to search for an "A" at the 5\'-end of the read and "TA" '
        "at 5'-overhang, use ``-m 'A,TA'``. Reads not containing the given motif are filtered out. ",
    ),
    "genome_2bit": typer.Option(
        "-g",
        "--genome-2bit",
        metavar=".2bit",
        default=None,
        metavar="file.2bit",
        help="If ``--motifFilter`` is provided, please also provide the genome sequence in 2bit format.",
    ),
    "gc_content_filter": typer.Option(
        "-gc",
        "--gc-content-filter",
        default=None,
        parser=parse_gc_content_filter,
        metavar="str",
        help="Check whether the GC content of the read falls within the provided range. Input format "
        "must be '<low>,<high>' , where <low> is the lower bound and <high> is the upper bound in "
        "the fraction of GC (eg. '0.1,0.9' ). If the GC content of the reads fall outside the range, "
        "they are filtered out. ",
    ),
    "min_aligned_fraction": typer.Option(
        "--min-aligned-fraction",
        default=None,
        metavar="float",
        help="Minimum fraction of the reads which should be aligned to be counted. This includes "
        "mismatches tolerated by the aligners, but excludes InDels/Clippings. (Default: %(default)s)",
    ),
}


READ_OPTS = {
    "min_mapping_quality": typer.Option(
        "-mq",
        "--min-mapping-quality",
        default=None,
        metavar="int",
        help="If set, only reads that have a mapping quality score of at least this are considered.",
    ),
    "sam_flag_include": typer.Option(
        "--sam-flag-include",
        default=None,
        metavar="int",
        help="Include reads based on SAM flag. For example, to get only reads that are the first "
        "mate, use a flag of 64. This is useful to count properly paired reads only once, as "
        "otherwise the second mate will be also considered for the coverage. This argument can be "
        "used more than once in a command. (Default: %(default)s)",
    ),
    "sam_flag_exclude": typer.Option(
        "--sam-flag-exclude",
        default=None,
        metavar="int",
        help="Exclude reads based on the SAM flag. For example, to get only reads that map to the "
        "forward strand, use ``--samFlagExclude 16``, where 16 is the SAM flag for reads that map "
        "to the reverse strand. This argument can be used more than once in a command. "
        "(Default: %(default)s)",
    ),
    "min_fragment_length": typer.Option(
        "--min-fragment-length",
        default=0,
        metavar="int",
        help="The minimum fragment length needed for fragment/read inclusion. This option is useful in "
        "ATACseq experiments, for filtering mono- or di-nucleosome fragments. (Default: %(default)s)",
    ),
    "max_fragment_length": typer.Option(
        "--max-fragment-length",
        default=None,
        metavar="int",
        help="The maximum fragment length accepted for fragment/read inclusion. (Default: %(default)s)",
    ),
    "filter_rna_strand": typer.Option(
        "--filter-rna-strand",
        metavar="str",
        default=FilterRNAStrand.none,
        help="Selects RNA-seq reads (single-end or paired-end) originating from genes on the given "
        "strand. This option assumes a standard dUTP-based library preparation (that is, "
        "``--filterRNAstrand forward`` keeps minus-strand reads, which originally came from genes"
        "on the forward strand using a dUTP-based method). Consider using "
        "``--samExcludeFlag`` instead for filtering by strand in other contexts.",
    ),
    "extend_reads": typer.Option(
        "-e",
        "--extendReads",
        default=None,
        metavar="int",
        help="This parameter allows the extension of reads to fragment size. If set, each read is "
        "extended, without exception.\n\n"
        "*NOTE*: This feature is generally NOT recommended for spliced-read data, such as "
        "RNA-seq, as it would extend reads over skipped regions.\n\n"
        "*Single-end*: Requires a user specified value for the final fragment length. Reads "
        "that already exceed this fragment length will not be extended.\n\n"
        "*Paired-end*: Reads with mates are always extended to match the fragment size defined "
        "by the two read mates. Unmated reads, mate reads that map too far apart (>4x fragment "
        "length) or even map to different chromosomes are treated like single-end reads. The "
        "input of a fragment length value is optional. If no value is specified, it is "
        "estimated from the data (mean of the fragment size of all mate reads).",
    ),
    "center_reads": typer.Option(
        "--center-reads",
        default=False,
        help="By adding this option, reads are centered with respect to the fragment length. For "
        "paired-end data, the read is centered at the fragment length defined by the two ends of"
        "the fragment. For single-end data, the given fragment length is used. This option is "
        "useful to get sharper signal around enriched regions.",
    )
}


### Enums for parameters with choices
class PlotFileFormat(str, Enum):
    png = "png"
    jpg = "jpg"
    svg = "svg"
    pdf = "pdf"


class DuplicateFilter(Enum):
    none = None
    start_bc = "start_bc"
    start__bc_umi = "start_umi"
    start_end_bc = "start_end_bc"
    start_end_bc_umi = "start_end_bc_umi"

class FilterRNAStrand(Enum):
    none = None
    forward = "forward"
    reverse = "reverse"


### Callbacks and parsers for options
def version_option() -> None:
    typer.echo(f"sincei {internal.version()}")
    raise typer.Exit()


def help_option(ctx: typer.Context) -> None:
    typer.echo(ctx.get_help())
    raise typer.Exit()


def parse_region(ctx: typer.Context, param: typer.Parameter, string: str | None) -> str | None:
    if string is None:
        return None

    # Remove whitespaces using split,join trick
    region = "".join(string.split())
    if region == "":
        return None

    # Remove undesired characters that may be present and replace "-" by ":"
    region = region.translate({ord(character): None for character in ",;|!{}()"}).replace("-", ":")

    if len(region) == 0:
        raise typer.BadParameter(f"{string} is not a valid region", param=param, ctx=ctx)

    return region


def parse_processors(ctx: typer.Context, param: typer.Parameter, value: str) -> int:
    # Use os.sched_getaffinity if available (Linux)
    # Fallback to os.cpu_count()
    try:
        available_processors = len(os.sched_getaffinity(0))
    except AttributeError:
        available_processors = os.cpu_count() or 1

    if value == "max/2":
        return int(available_processors * 0.5)
    if value == "max":
        return available_processors

    try:
        number_of_processors = int(value)
    except ValueError as exc:
        raise typer.BadParameter(f"{value} is not a valid number of processors", param=param, ctx=ctx) from exc

    if number_of_processors > available_processors:
        number_of_processors = available_processors

    return number_of_processors


def parse_motif_filter(ctx: typer.Context, param: typer.Parameter, value: str) -> tuple[str, str]:
    motifs = value.strip(" ").split(",")
    if len(motifs) != 2:
        raise typer.BadParameter(
            f"{value} is not a valid motif filter. Please provide two motifs separated by a comma, for example 'A,TA'.",
            param=param,
            ctx=ctx,
        )
    return motifs[0], motifs[1]


def parse_gc_content_filter(ctx: typer.Context, param: typer.Parameter, value: str) -> tuple[float, float]:
    gc = value.strip(" ").split(",")
    return float(gc[0]), float(gc[1])


# Typer does not support whitespace-separated multi-value options, so we preprocess the sysargv
# to add the option flag before each value. This is a workaround to allow for more user-friendly
# CLI commands where options can be followed by multiple values without needing to repeat the option
# flag for each value. The downside is that options should always be after arguments in the CLI
# command for this to work correctly.
def preprocess_args():
    """
    We preprocess the sysargv so that:
    - python3 app.py some_command --filters filter1 filter2 filter3 --environments env1 env2 env3

    becomes:
    - python3 app.py some_command --filters filter1 --filters filter2 --filters filter3 --environments env1 --environments env2 --environments env3

    //!\\ DOWNSIDE: options should always be after arguments in the CLI command //!\\
    
    We do not care about this sincei we only use options.
    """

    # Get main cmd
    final_cmd = []
    for idx, arg in enumerate(sys.argv):
        if any(arg.startswith(_) for _ in ["-", "--"]):
            break
        else:
            final_cmd.append(arg)

    # Get options and their values
    for idx, arg in enumerate(sys.argv):
        if any(arg.startswith(_) for _ in ["-", "--"]):
            opt_values = []
            for value in sys.argv[idx + 1 :]:
                if any(value.startswith(_) for _ in ["-", "--"]):
                    break
                else:
                    opt_values.append(value)

            if len(opt_values) >= 1:
                [final_cmd.extend([arg, opt_value]) for opt_value in opt_values]
            else:
                final_cmd.append(arg)

    # Replace by reformatted command
    sys.argv = final_cmd


def validateInputs(
    args: argparse.Namespace, process_barcodes: bool = True
) -> argparse.Namespace | tuple[argparse.Namespace, list[str]]:
    """
    Ensure that right input is provided from argparse.
    """

    ## Labels
    if args.groupTag:
        # in case of --groupTag, use args.labels as groups
        if len(args.bamfiles) > 1:
            logging.info("Only a single BAM file is allowed when --groupTag is specified.")
            exit(0)
        if not args.labels:
            logging.info("Please indicate the sample groups to be processed from the BAM file with --labels")
            exit(0)
    else:
        if args.labels and len(args.bamfiles) != len(args.labels):
            logging.info(
                "The number of labels does not match the number of bam files. "
                "This is only allowed if a single BAM file is provided and --groupTag is specified."
            )
            exit(0)
        if not args.labels:
            args.labels = smartLabels(args.bamfiles)

    ## Motif and GC filter
    if args.motifFilter:
        if not args.genome2bit:
            logging.info("MotifFilter asked but genome (2bit) file not provided.")
            sys.exit(1)
        else:
            args.motifFilter = [x.strip(" ").split(",") for x in args.motifFilter]

    if args.GCcontentFilter:
        gc = args.GCcontentFilter.strip(" ").split(",")
        args.GCcontentFilter = [float(x) for x in gc]

    ## Barcodes + new labels (if required)
    if process_barcodes:
        with open(args.barcodes, "r") as f:
            barcodes = f.read().splitlines()
        f.close()
        args.barcodes = barcodes
        newlabels = ["{}::{}".format(a, b) for a in args.labels for b in barcodes]
        return args, newlabels
    else:
        return args