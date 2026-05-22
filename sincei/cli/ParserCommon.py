# -*- coding: utf-8 -*-
from __future__ import annotations

import argparse
import os
import sys
import logging

from sincei import _sincei as internal


def inputOutputOptions(opts=None, requiredOpts=[], suppress_args=None) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group("Input/Output options")

    ## inputs
    if "h5adfile" in opts:
        group.add_argument(
            "-i",
            "--input",
            metavar=".h5ad",
            help="sincei-generated input file in .h5ad format.",
            required=True,
        )
    elif "h5adfiles" in opts:
        group.add_argument(
            "-i",
            "--input",
            metavar=".h5ad",
            help="List of sincei-generated input .h5ad files separated by spaces.",
            nargs="+",
            required=True,
        )

    if "h5mufile" in opts:
        group.add_argument(
            "-i",
            "--input",
            metavar=".h5mu",
            help="sincei-generated input file in .h5mu format.",
            required=True,
        )

    if "bamfiles" in opts:
        group.add_argument(
            "-b",
            "--bamfiles",
            metavar=".bam",
            help="List of indexed BAM files separated by spaces.",
            nargs="+",
            required=True,
        )
    if "bamfile" in opts:
        group.add_argument(
            "-b",
            "--bamfile",
            metavar=".bam",
            help="Indexed BAM file.",
            required=True,
        )

    if "whitelist" in opts:
        group.add_argument(
            "-w",
            "--whitelist",
            metavar=".txt",
            help="A single-column file containing the whitelist of barcodes to be used.",
            default=None,
            required=True if "whitelist" in requiredOpts else False,
        )
    if "barcodes" in opts:
        group.add_argument(
            "-bc",
            "--barcodes",
            metavar=".txt",
            help="A single-column file containing the cell barcode whitelist, one barcode per line.",
            required=True if "barcodes" in requiredOpts else False,
        )

    if "groupInfo" in opts:
        group.add_argument(
            "-i",
            "--groupInfo",
            metavar=".tsv",
            help="A 4-column tsv file with cell grouping information in the "
            "format: `sample::barcode, UMAP1, UMAP2, group` (like the output from scClusterCells) "
            "or 3-column tsv file with format: `sample, barcode, group`. "
            "Coverages will be computed per group.",
            required=True,
        )

    if "BED" in opts:
        group.add_argument(
            "--BED",
            metavar=(".bed", ".gtf"),
            help=show_or_hide(
                "BED/GTF files to limit the coverage analysis to the regions in them.",
                "BED",
                suppress_args,
            ),
            nargs="+",
            required=True if "BED" in requiredOpts else False,
        )

    ## outputs
    if "outFilePrefix" in opts:
        group.add_argument(
            "-o",
            "--outFilePrefix",
            metavar="STR",
            help="Prefix for output file names.",
            type=str,
            required=True if "outFilePrefix" in requiredOpts else False,
        )

    if "outFile" in opts:
        group.add_argument(
            "-o",
            "--outFile",
            metavar="STR",
            help="The file to write results to. For `scFilterStats`, `scFilterBarcodes` "
            "and `scJSD`, the output file is a .tsv file. For other tools, the output file is "
            "an updated .h5ad or .h5mu file with the result of the requested operation.",
            type=str,
            required=True if "outFile" in requiredOpts else False,
        )

    if "labels" in opts:
        group.add_argument(
            "-l",
            "--labels",
            metavar="sample",
            help="User defined labels instead of default labels from file names. "
            "Multiple labels have to be separated by a space, e.g. "
            "``--labels sample1 sample2 sample3``.",
            nargs="+",
        )

        group.add_argument(
            "--smartLabels",
            action="store_true",
            help="Instead of manually specifying labels for the input BAM files, use the "
            "file name after removing the path and extension.",
        )

    if "region" in opts:
        group.add_argument(
            "-r",
            "--region",
            metavar="CHR:START:END",
            help="Region of the genome to limit the operation to - this is useful when "
            "testing parameters to reduce the computing time. The format is "
            "chr:start:end, for example ``--region chr10`` or ``--region chr10:456700:891000``.",
            required=False,
            type=genomicRegion,
        )

    return parser


def otherOptions() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group("Other options")

    group.add_argument(
        "-p",
        "--numberOfProcessors",
        metavar="INT",
        help='Number of processors to use. Type "max/2" to use half the maximum number of '
        'processors or "max" to use all available processors. (Default: "max")',
        type=numberOfProcessors,
        default=numberOfProcessors("max"),
        required=False,
    )

    group.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this message and exit",
    )

    group.add_argument(
        "-v",
        "--verbose",
        help="Set to see processing messages.",
        action="store_true",
    )

    group.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"%(prog)s {internal.version()}",
        help="Print the program version and exit.",
    )

    return parser


def gtf_options() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group("GTF/BED12 options")

    group.add_argument(
        "--metagene",
        help="When either a BED12 or GTF file are used to provide regions, "
        "perform the computation on the merged exons, rather than using the "
        "genomic interval defined by the 5-prime and 3-prime most transcript "
        "bound (i.e., columns 2 and 3 of a BED file). If a BED3 or BED6 file is "
        "used as input, then columns 2 and 3 are used as an exon. "
        "(Default: %(default)s)",
        action="store_true",
        dest="keepExons",
    )

    group.add_argument(
        "--transcriptID",
        metavar="STR",
        help="When a GTF file is used to provide regions, only  entries with "
        "this value as their feature (column 3) will be processed as "
        "transcripts. (Default: %(default)s)",
        default="transcript",
    )

    group.add_argument(
        "--exonID",
        metavar="STR",
        help="When a GTF file is used to provide regions, only entries with "
        "this value as their feature (column 3) will be processed as exons. "
        "CDS would be another common value for this. (Default: %(default)s)",
        default="exon",
    )

    group.add_argument(
        "--transcript_id_designator",
        metavar="STR",
        help="Each region has an ID (e.g., ACTB) assigned to it, which for BED "
        "files is either column 4 (if it exists) or the interval bounds. For "
        "GTF files this is instead stored in the last column as a key:value pair "
        "(e.g., as 'transcript_id \"ACTB\"', for a key of transcript_id and a "
        "value of ACTB). In some cases it can be convenient to use a different "
        "identifier. To do so, set this to the desired key. (Default: %(default)s)",
        default="transcript_id",
    )

    return parser


def plotOptions() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group("Plot options")

    group.add_argument(
        "--plotWidth",
        metavar="INT",
        default=10,
        type=float,
        help="Output plot width (in cm). (Default: %(default)s)",
    )

    group.add_argument(
        "--plotHeight",
        metavar="INT",
        default=10,
        type=float,
        help="Output plot height (in cm). (Default: %(default)s)",
    )

    group.add_argument(
        "--plotFileFormat",
        help="Image format type. If given, this option "
        "overrides the image format based on the plotFile "
        "ending. (Default: %(default)s)",
        type=str,
        choices=["png", "pdf", "svg", "eps"],
        default="png",
    )

    return parser


## Tools: scBulkCoverage, scCountReads, scFilterStats, scJSD
def bamOptions(suppress_args=None, default_opts=None) -> argparse.ArgumentParser:
    """
    Arguments for tools that operate on BAM files.
    """

    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group("BAM processing options")

    group.add_argument(
        "-ct",
        "--cellTag",
        metavar="STR",
        help="Name of the BAM tag from which to extract barcodes. (Default: %(default)s)",
        type=str,
        default="BC",
    )

    group.add_argument(
        "-gt",
        "--groupTag",
        metavar="STR",
        help=show_or_hide(
            "In case of a groupped BAM file, such as the one containing Read Group (``RG``) or Sample (``SM``) tag,"
            "it is possible to process group the reads using the provided ``--groupTag`` argument. NOTE: In case of such input, "
            "please ensure that the ``--labels`` argument indicates the expected group labels contained in the BAM files. "
            "The ``--groupTag`` along with the ``--cellTag`` is then used to identify unique samples (cells) from the input.",
            "groupTag",
            suppress_args,
        ),
        type=str,
        default=None,
    )

    group.add_argument(
        "-l",
        "--labels",
        metavar="sample",
        help=show_or_hide(
            "User defined labels instead of default labels from file names. "
            "Multiple labels have to be separated by a space, e.g. "
            "``--labels sample1 sample2 sample3``.",
            "labels",
            suppress_args,
        ),
        nargs="+",
    )

    group.add_argument(
        "--smartLabels",
        action="store_true",
        help=show_or_hide(
            "Instead of manually specifying labels for the input BAM files, use "
            "the file name after removing the path and extension.",
            "smartLabels",
            suppress_args,
        ),
    )

    group.add_argument(
        "-r",
        "--region",
        metavar="CHR:START:END",
        help=show_or_hide(
            "Region of the genome to limit the operation to - this is useful when "
            "testing parameters to reduce the computing time. The format is "
            "chr:start:end, for example ``--region chr10`` or ``--region chr10:456700:891000``.",
            "region",
            suppress_args,
        ),
        required=False,
        type=genomicRegion,
    )

    group.add_argument(
        "-bl",
        "--blacklist",
        "--blackListFileName",
        metavar=".bed",
        help=show_or_hide(
            "A BED or GTF file containing regions that should be excluded from all analyses. "
            "Currently this works by rejecting genomic chunks that happen to overlap an entry. "
            "Consequently, for BAM files, if a read partially overlaps a blacklisted region or a "
            "fragment spans over it, then the read/fragment might still be considered. Please note "
            "that you should adjust the effective genome size, if relevant.",
            "blacklist",
            suppress_args,
        ),
        nargs="+",
        dest="blackListFileName",
        required=False,
    )

    group.add_argument(
        "--chrToSkip",
        metavar="CHR",
        help=show_or_hide(
            "List of chromosome names to skip from the analysis. "
            "Regions on these chromosomes will be excluded. "
            "Useful for skipping mitochondrial, X chromosome, or unplaced contigs. "
            "Multiple chromosomes can be specified, e.g. ``--chrToSkip chrM chrX``.",
            "chrToSkip",
            suppress_args,
        ),
        nargs="+",
        default=None,
        required=False,
    )

    group.add_argument(
        "-bs",
        "--binSize",
        metavar="INT",
        help=show_or_hide(
            "Size of the bins, in bases, to calculate coverage. (Default: %(default)s)",
            "binSize",
            suppress_args,
        ),
        type=int,
        default=get_default("binSize", default_opts),
    )

    group.add_argument(
        "--distanceBetweenBins",
        metavar="INT",
        help=show_or_hide(
            "The gap length, in bases, between bins for calculating coverage. "
            "Larger values can be used to sample the genome for input files with high coverage "
            "while smaller values are useful for lower coverage data. Note that if you specify a value that "
            "results in too few (<1000) reads sampled, the value will be decreased. (Default: %(default)s)",
            "distanceBetweenBins",
            suppress_args,
        ),
        default=get_default("distanceBetweenBins", default_opts),
        type=int,
    )

    return parser


## Tools: scBulkCoverage, scCountReads
def filterOptions(suppress_args=None) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group("Read Filtering Options")

    group.add_argument(
        "--duplicateFilter",
        help=show_or_hide(
            "How to filter for duplicates? Different combinations (using start/end/umi) are possible. "
            "Read start position and read barcode are always considered. Default (None) considers "
            "all reads as passing the filter. "
            "Note that in case of paired end data, both reads in the fragment are considered (and kept). "
            "So if you wish to keep only read1, combine this option with `--samFlagInclude`. ",
            "duplicateFilter",
            suppress_args,
        ),
        type=str,
        choices=["start_bc", "start_bc_umi", "start_end_bc", "start_end_bc_umi"],
        default=None,
    )

    group.add_argument(
        "-m",
        "--motifFilter",
        metavar="STR",
        help=show_or_hide(
            "Check whether a given motif is present in the read and the corresponding reference genome. "
            "This option checks for the motif at the 5-end of the read and at the 5-overhang in the genome, "
            "which is useful in identifying reads properly cut by a restriction-enzyme or MNAse. "
            'For example, if you want to search for an "A" at the 5\'-end of the read and "TA" at 5\'-overhang, '
            "use ``-m 'A,TA'``. Reads not containing the given motif are filtered out. ",
            "motifFilter",
            suppress_args,
        ),
        type=str,
        nargs="+",
        default=None,
    )

    group.add_argument(
        "-g",
        "--genome2bit",
        metavar=".2bit",
        help=show_or_hide(
            "If ``--motifFilter`` is provided, please also provide the genome sequence in 2bit format.",
            "genome2bit",
            suppress_args,
        ),
        type=str,
        default=None,
    )

    group.add_argument(
        "-gc",
        "--GCcontentFilter",
        metavar="STR",
        help=show_or_hide(
            "Check whether the GC content of the read falls within the provided range "
            "Input format must be '<low>,<high>' , where <low> is the lower bound and <high> is the "
            "upper bound in the fraction of GC (eg. '0.1,0.9' ). "
            "If the GC content of the reads fall outside the range, they are filtered out. ",
            "GCcontentFilter",
            suppress_args,
        ),
        type=str,
        default=None,
    )

    group.add_argument(
        "--minAlignedFraction",
        metavar="FLOAT",
        help=show_or_hide(
            "Minimum fraction of the reads which should be aligned to be counted. This includes "
            "mismatches tolerated by the aligners, but excludes InDels/Clippings. (Default: %(default)s)",
            "minAlignedFraction",
            suppress_args,
        ),
        default=None,
        type=float,
        required=False,
    )

    return parser


## Tools: scBulkCoverage, scCountReads
def readOptions(suppress_args=None) -> argparse.ArgumentParser:
    """Common arguments related to BAM files and the interpretation
    of the read coverage
    """
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group("Read Processing Options")

    group.add_argument(
        "--minMappingQuality",
        metavar="INT",
        help=show_or_hide(
            "If set, only reads that have a mapping quality score of at least this are considered.",
            "minMappingQuality",
            suppress_args,
        ),
        type=int,
    )

    group.add_argument(
        "--samFlagInclude",
        metavar="INT",
        help=show_or_hide(
            "Include reads based on SAM flag. For example, to get only reads that are the first "
            "mate, use a flag of 64. This is useful to count properly paired reads only once, as "
            "otherwise the second mate will be also considered for the coverage. This argument can "
            "be used more than once in a command. (Default: %(default)s)",
            "samFlagInclude",
            suppress_args,
        ),
        default=None,
        type=int,
        required=False,
    )

    group.add_argument(
        "--samFlagExclude",
        metavar="INT",
        help=show_or_hide(
            "Exclude reads based on the SAM flag. For example, to get only reads that map to the "
            "forward strand, use ``--samFlagExclude 16``, where 16 is the SAM flag for reads that "
            "map to the reverse strand. This argument can be used more than once in a command. "
            "(Default: %(default)s)",
            "samFlagExclude",
            suppress_args,
        ),
        default=None,
        type=int,
        required=False,
    )

    group.add_argument(
        "--minFragmentLength",
        metavar="INT",
        help=show_or_hide(
            "The minimum fragment length needed for read/pair inclusion. This option is for "
            "useful in ATACseq experiments, for filtering mono- or di-nucleosome fragments. "
            "(Default: %(default)s)",
            "minFragmentLength",
            suppress_args,
        ),
        default=0,
        type=int,
        required=False,
    )

    group.add_argument(
        "--maxFragmentLength",
        metavar="INT",
        help=show_or_hide(
            "The maximum fragment length accepted for read/pair inclusion. (Default: %(default)s)",
            "maxFragmentLength",
            suppress_args,
        ),
        default=0,
        type=int,
        required=False,
    )

    group.add_argument(
        "--filterRNAstrand",
        help=show_or_hide(
            "Selects RNA-seq reads (single-end or paired-end) originating from genes on the given "
            "strand. This option assumes a standard dUTP-based library preparation (that is, "
            "``--filterRNAstrand=forward`` keeps minus-strand reads, which originally came from "
            "genes on the forward strand using a dUTP-based method). Consider using "
            "``--samExcludeFlag`` instead for filtering by strand in other contexts.",
            "filterRNAstrand",
            suppress_args,
        ),
        choices=["forward", "reverse"],
        default=None,
    )

    group.add_argument(
        "-e",
        "--extendReads",
        metavar="INT",
        help=show_or_hide(
            "This parameter allows the extension of reads to fragment size. If set, each read is "
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
            "extendReads",
            suppress_args,
        ),
        type=int,
        nargs="?",
        const=True,
        default=False,
    )

    group.add_argument(
        "--centerReads",
        help=show_or_hide(
            "By adding this option, reads are centered with respect to the fragment length. For "
            "paired-end data, the read is centered at the fragment length defined by the two ends "
            "of the fragment. For single-end data, the given fragment length is used. This option "
            "is useful to get sharper signal around enriched regions.",
            "centerReads",
            suppress_args,
        ),
        action="store_true",
    )

    return parser


## helper functions
def show_or_hide(help: str, name: str, arglist: list[str] | None) -> str:
    if arglist and name in arglist:
        return argparse.SUPPRESS
    else:
        return help


def get_default(name: str, argdict: dict[str, object] | None) -> object | None:
    if argdict and name in argdict.keys():
        return argdict[name]
    else:
        return None


def process_args(args: argparse.Namespace) -> argparse.Namespace:
    if not args.labels:
        if args.smartLabels:
            args.labels = smartLabels(args.bamfiles)
        else:
            args.labels = args.bamfiles

    if len(args.bamfiles) != len(args.labels):
        sys.exit("The number of labels does not match the number of BAM files.")

    return args


def genomicRegion(string: str) -> str | None:
    # remove whitespaces using split,join trick
    region = "".join(string.split())
    if region == "":
        return None
    # remove undesired characters that may be present and
    # replace - by :
    # N.B., the syntax for translate() differs between python 2 and 3
    try:
        region = region.translate(None, ",;|!{}()").replace("-", ":")
    except:
        region = region.translate({ord(i): None for i in ",;|!{}()"})
    if len(region) == 0:
        raise argparse.ArgumentTypeError("{} is not a valid region".format(string))
    return region


def numberOfProcessors(string: str) -> int:
    try:
        # Use os.sched_getaffinity if available (Linux)
        availProc = len(os.sched_getaffinity(0))
    except AttributeError:
        # Fallback to os.cpu_count()
        availProc = os.cpu_count()

    if string == "max/2":
        numberOfProcessors = int(availProc * 0.5)
    elif string == "max":  # default case
        # by default, use all available processors
        numberOfProcessors = availProc
    else:
        try:
            numberOfProcessors = int(string)
        except ValueError:
            raise argparse.ArgumentTypeError(f"{string} is not a valid number of processors")

        except Exception as e:
            raise argparse.ArgumentTypeError(
                f"The given value {string} is not valid. Error message: {e}\n"
                f"The number of available processors in your computer is {availProc}."
            )

        if numberOfProcessors > availProc:
            numberOfProcessors = availProc

    return numberOfProcessors


def smartLabel(label: str) -> str:
    """
    Remove the path name and the last extension from the file name.
    /pth/to/file.name.bam -> file.name
    """
    lab = os.path.splitext(os.path.basename(label))[0]
    if lab == "":
        # Maybe we have a dot file?
        lab = os.path.basename(label)
    return lab


def smartLabels(labels: list[str]) -> list[str]:
    smrt = [smartLabel(x) for x in labels]
    if len(smrt) != len(set(smrt)):
        print(
            "Labels inferred from file names are not unique. "
            "Please be aware that in case of overlapping barcodes the counts will be merged."
        )
    return smrt


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
