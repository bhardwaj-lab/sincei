from __future__ import annotations

import typer

from .shared import log_parameters, normalize_processors, normalize_region, version_option

DESCRIPTION = """
Get pseudo-bulk coverage per group using a cell->group mapping (output of scClusterCells).
"""


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help=DESCRIPTION,
    context_settings={"help_option_names": ["-h", "--help"]},
)


@version_option("scBulkCoverage")
@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    bamfiles: list[str] = typer.Option(..., "-b", "--bamfiles", help="List of indexed BAM files separated by spaces."),
    groupInfo: str = typer.Option(
        ...,
        "-i",
        "--groupInfo",
        help="A 4-column tsv file with cell grouping information in the format: `sample::barcode, UMAP1, UMAP2, group` (like the output from scClusterCells) or 3-column tsv file with format: `sample, barcode, group`. Coverages will be computed per group.",
    ),
    outFilePrefix: str = typer.Option(..., "-o", "--outFilePrefix", help="Prefix for output file names."),
    outFileFormat: str = typer.Option(
        "bigwig", "-of", "--outFileFormat", help='Output file type. Either "bigwig" or "bedgraph".'
    ),
    normalizeUsing: str = typer.Option(
        "CPM",
        "-n",
        "--normalizeUsing",
        help='How to normalize the pseudo-bulk counts. Options are "CPM": normalized each bin to the counts per million mapped reads in that group.\n"Frequency": binarize the coverage per bin and normalize to the total no. of cells per group. \n"Mean": get mean signal per bin across cells in each group.\n"None": simply return the sum of coverage per group.',
    ),
    ignoreForNormalization: list[str] = typer.Option(
        None,
        "-ig",
        "--ignoreForNormalization",
        metavar="CHR",
        help="Chromosomes to skip while calculating normalization factors",
    ),
    normalizeByReference: str = typer.Option(
        None,
        "-nr",
        "--normalizeByReference",
        help="NOT IMPLEMENTED: Normalize each group of cells by a reference group (which must be present in the --groupinfo file)Note that the --normalizeUsing method is applied beforehand.",
    ),
    scaleFactor: float = typer.Option(
        1.0,
        "--scaleFactor",
        help="The computed scaling factor (or 1, if not applicable) will be multiplied by this. (Default: %(default)s)",
    ),
    MNase: bool = typer.Option(
        False,
        "--MNase",
        help="Determine nucleosome positions from MNase-seq/CUTnRUN data. Only 3 nucleotides at the center of each fragment are counted. The fragment ends are defined by the two mate reads. Only fragment lengthsbetween 130 - 200 bp are considered to avoid dinucleosomes or other artifacts. By default, any fragments smaller or larger than this are ignored. To over-ride this, use the --minFragmentLength and --maxFragmentLength options, which will default to 130 and 200 if not otherwise specified in the presence of --MNase. *NOTE*: Requires paired-end data. A bin size of 1 is recommended.",
    ),
    Offset: list[int] = typer.Option(
        None,
        "--Offset",
        help="Uses this offset inside of each read as the signal. This is useful in cases like RiboSeq or GROseq, where the signal is 12, 15 or 0 bases past the start of the read. This can be paired with the --filterRNAstrand option. Note that negative values indicate offsets from the end of each read. A value of 1 indicates the first base of the alignment (taking alignment orientation into account). Likewise, a value of -1 is the last base of the alignment. An offset of 0 is not permitted. If two values are specified, then they will be used to specify a range of positions. Note that specifying something like --Offset 5 -1 will result in the 5th through last position being used, which is equivalent to trimming 4 bases from the 5-prime end of alignments. Note that if you specify --centerReads, the centering will be performed before the offset.",
    ),
    cellTag: str = typer.Option(
        "BC", "-ct", "--cellTag", help="Name of the BAM tag from which to extract barcodes. (Default: %(default)s)"
    ),
    groupTag: str = typer.Option(
        None,
        "-gt",
        "--groupTag",
        help="In case of a groupped BAM file, such as the one containing Read Group (``RG``) or Sample (``SM``) tag,it is possible to process group the reads using the provided ``--groupTag`` argument. NOTE: In case of such input, please ensure that the ``--labels`` argument indicates the expected group labels contained in the BAM files. The ``--groupTag`` along with the ``--cellTag`` is then used to identify unique samples (cells) from the input.",
    ),
    labels: list[str] = typer.Option(
        None,
        "-l",
        "--labels",
        metavar="sample",
        help="User defined labels instead of default labels from file names. Multiple labels have to be separated by a space, e.g. ``--labels sample1 sample2 sample3``.",
    ),
    smartLabels: bool = typer.Option(
        False,
        "--smartLabels",
        help="Instead of manually specifying labels for the input BAM files, use the file name after removing the path and extension.",
    ),
    region: str = typer.Option(
        None,
        "-r",
        "--region",
        callback=normalize_region,
        help="Region of the genome to limit the operation to - this is useful when testing parameters to reduce the computing time. The format is chr:start:end, for example ``--region chr10`` or ``--region chr10:456700:891000``.",
    ),
    blackListFileName: list[str] = typer.Option(
        None,
        "-bl",
        "--blacklist",
        "--blackListFileName",
        metavar=".bed",
        help="A BED or GTF file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry. Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered. Please note that you should adjust the effective genome size, if relevant.",
    ),
    chrToSkip: list[str] = typer.Option(
        None,
        "--chrToSkip",
        metavar="CHR",
        help="List of chromosome names to skip from the analysis. Regions on these chromosomes will be excluded. Useful for skipping mitochondrial, X chromosome, or unplaced contigs. Multiple chromosomes can be specified, e.g. ``--chrToSkip chrM chrX``.",
    ),
    binSize: int = typer.Option(
        100, "-bs", "--binSize", help="Size of the bins, in bases, to calculate coverage. (Default: %(default)s)"
    ),
    distanceBetweenBins: int = typer.Option(
        None,
        "--distanceBetweenBins",
        help="The gap length, in bases, between bins for calculating coverage. Larger values can be used to sample the genome for input files with high coverage",
    ),
    duplicateFilter: str = typer.Option(
        None,
        "--duplicateFilter",
        help="How to filter for duplicates? Different combinations (using start/end/umi) are possible. Read start position and read barcode are always considered. Default (None) considers all reads as passing the filter. Note that in case of paired end data, both reads in the fragment are considered (and kept). So if you wish to keep only read1, combine this option with `--samFlagInclude`. ",
    ),
    motifFilter: list[str] = typer.Option(
        None,
        "-m",
        "--motifFilter",
        metavar="STR",
        help="Check whether a given motif is present in the read and the corresponding reference genome. This option checks for the motif at the 5-end of the read and at the 5-overhang in the genome, which is useful in identifying reads properly cut by a restriction-enzyme or MNAse. For example, if you want to search for an \"A\" at the 5'-end of the read and \"TA\" at 5'-overhang, use ``-m 'A,TA'``. Reads not containing the given motif are filtered out. ",
    ),
    genome2bit: str = typer.Option(
        None,
        "-g",
        "--genome2bit",
        metavar=".2bit",
        help="If ``--motifFilter`` is provided, please also provide the genome sequence in 2bit format.",
    ),
    GCcontentFilter: str = typer.Option(
        None,
        "-gc",
        "--GCcontentFilter",
        metavar="STR",
        help="Check whether the GC content of the read falls within the provided range Input format must be '<low>,<high>' , where <low> is the lower bound and <high> is the upper bound in the fraction of GC (eg. '0.1,0.9' ). If the GC content of the reads fall outside the range, they are filtered out. ",
    ),
    minAlignedFraction: float = typer.Option(
        None,
        "--minAlignedFraction",
        metavar="FLOAT",
        help="Minimum fraction of the reads which should be aligned to be counted. This includes mismatches tolerated by the aligners, but excludes InDels/Clippings. (Default: %(default)s)",
    ),
    minMappingQuality: int = typer.Option(
        None,
        "--minMappingQuality",
        metavar="INT",
        help="If set, only reads that have a mapping quality score of at least this are considered.",
    ),
    samFlagInclude: int = typer.Option(
        None,
        "--samFlagInclude",
        metavar="INT",
        help="Include reads based on SAM flag. For example, to get only reads that are the first mate, use a flag of 64. This is useful to count properly paired reads only once, as otherwise the second mate will be also considered for the coverage. This argument can be used more than once in a command. (Default: %(default)s)",
    ),
    samFlagExclude: int = typer.Option(
        None,
        "--samFlagExclude",
        metavar="INT",
        help="Exclude reads based on the SAM flag. For example, to get only reads that map to the forward strand, use ``--samFlagExclude 16``, where 16 is the SAM flag for reads that map to the reverse strand. This argument can be used more than once in a command. (Default: %(default)s)",
    ),
    minFragmentLength: int = typer.Option(
        0,
        "--minFragmentLength",
        metavar="INT",
        help="The minimum fragment length needed for read/pair inclusion. This option is for useful in ATACseq experiments, for filtering mono- or di-nucleosome fragments. (Default: %(default)s)",
    ),
    maxFragmentLength: int = typer.Option(
        0,
        "--maxFragmentLength",
        metavar="INT",
        help="The maximum fragment length accepted for read/pair inclusion. (Default: %(default)s)",
    ),
    filterRNAstrand: str = typer.Option(
        None,
        "--filterRNAstrand",
        help="Selects RNA-seq reads (single-end or paired-end) originating from genes on the given strand. This option assumes a standard dUTP-based library preparation (that is, ``--filterRNAstrand=forward`` keeps minus-strand reads, which originally came from genes on the forward strand using a dUTP-based method). Consider using ``--samExcludeFlag`` instead for filtering by strand in other contexts.",
    ),
    extendReads: int = typer.Option(
        False,
        "-e",
        "--extendReads",
        metavar="INT",
        flag_value=True,
        help="This parameter allows the extension of reads to fragment size. If set, each read is extended, without exception.\n\n*NOTE*: This feature is generally NOT recommended for spliced-read data, such as RNA-seq, as it would extend reads over skipped regions.\n\n*Single-end*: Requires a user specified value for the final fragment length. Reads that already exceed this fragment length will not be extended.\n\n*Paired-end*: Reads with mates are always extended to match the fragment size defined by the two read mates. Unmated reads, mate reads that map too far apart (>4x fragment length) or even map to different chromosomes are treated like single-end reads. The input of a fragment length value is optional. If no value is specified, it is estimated from the data (mean of the fragment size of all mate reads).",
    ),
    centerReads: bool = typer.Option(
        False,
        "--centerReads",
        help="By adding this option, reads are centered with respect to the fragment length. For paired-end data, the read is centered at the fragment length defined by the two ends of the fragment. For single-end data, the given fragment length is used. This option is useful to get sharper signal around enriched regions.",
    ),
    numberOfProcessors: str = typer.Option(
        "max",
        "-p",
        "--numberOfProcessors",
        callback=normalize_processors,
        help='Number of processors to use. Type "max/2" to use half the maximum number of processors or "max" to use all available processors. (Default: "max")',
    ),
    verbose: bool = typer.Option(False, "-v", "--verbose", help="Set to see processing messages."),
) -> int:
    if ctx.invoked_subcommand is None:
        log_parameters(
            bamfiles=bamfiles,
            groupInfo=groupInfo,
            outFilePrefix=outFilePrefix,
            outFileFormat=outFileFormat,
            normalizeUsing=normalizeUsing,
            ignoreForNormalization=ignoreForNormalization,
            normalizeByReference=normalizeByReference,
            scaleFactor=scaleFactor,
            MNase=MNase,
            Offset=Offset,
            cellTag=cellTag,
            groupTag=groupTag,
            labels=labels,
            smartLabels=smartLabels,
            region=region,
            blackListFileName=blackListFileName,
            chrToSkip=chrToSkip,
            binSize=binSize,
            distanceBetweenBins=distanceBetweenBins,
            duplicateFilter=duplicateFilter,
            motifFilter=motifFilter,
            genome2bit=genome2bit,
            GCcontentFilter=GCcontentFilter,
            minAlignedFraction=minAlignedFraction,
            minMappingQuality=minMappingQuality,
            samFlagInclude=samFlagInclude,
            samFlagExclude=samFlagExclude,
            minFragmentLength=minFragmentLength,
            maxFragmentLength=maxFragmentLength,
            filterRNAstrand=filterRNAstrand,
            extendReads=extendReads,
            centerReads=centerReads,
            numberOfProcessors=numberOfProcessors,
            verbose=verbose,
        )
    return 0


if __name__ == "__main__":
    app()
