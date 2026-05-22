from __future__ import annotations

import typer

from .shared import log_parameters, normalize_processors, normalize_region, version_option

DESCRIPTION = """
Counts reads for each cell barcode on genomic bins or user-defined features.
"""


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help=DESCRIPTION,
    context_settings={"help_option_names": ["-h", "--help"]},
)

bins_app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Count reads in fixed-size bins.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
features_app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Count reads on provided features.",
    context_settings={"help_option_names": ["-h", "--help"]},
)

app.add_typer(bins_app, name="bins")
app.add_typer(features_app, name="features")


@app.callback(invoke_without_command=True)
def main(ctx: typer.Context) -> int:
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
        raise typer.Exit()
    return 0


def _count_reads(
    bamfiles: list[str],
    barcodes: str,
    outFile: str,
    GTF: str,
    labels: list[str],
    smartLabels: bool,
    region: str,
    blackListFileName: list[str],
    chrToSkip: list[str],
    binSize: int,
    distanceBetweenBins: int,
    cellTag: str,
    groupTag: str,
    minMappingQuality: int,
    samFlagInclude: int,
    samFlagExclude: int,
    minFragmentLength: int,
    maxFragmentLength: int,
    filterRNAstrand: str,
    extendReads: int,
    centerReads: bool,
    duplicateFilter: str,
    motifFilter: list[str],
    genome2bit: str,
    GCcontentFilter: str,
    minAlignedFraction: float,
    valueTag: str,
    genomeChunkSize: int,
    numberOfProcessors: str,
    verbose: bool,
) -> int:
    log_parameters(
        bamfiles=bamfiles,
        barcodes=barcodes,
        outFile=outFile,
        GTF=GTF,
        labels=labels,
        smartLabels=smartLabels,
        region=region,
        blackListFileName=blackListFileName,
        chrToSkip=chrToSkip,
        binSize=binSize,
        distanceBetweenBins=distanceBetweenBins,
        cellTag=cellTag,
        groupTag=groupTag,
        minMappingQuality=minMappingQuality,
        samFlagInclude=samFlagInclude,
        samFlagExclude=samFlagExclude,
        minFragmentLength=minFragmentLength,
        maxFragmentLength=maxFragmentLength,
        filterRNAstrand=filterRNAstrand,
        extendReads=extendReads,
        centerReads=centerReads,
        duplicateFilter=duplicateFilter,
        motifFilter=motifFilter,
        genome2bit=genome2bit,
        GCcontentFilter=GCcontentFilter,
        minAlignedFraction=minAlignedFraction,
        valueTag=valueTag,
        genomeChunkSize=genomeChunkSize,
        numberOfProcessors=numberOfProcessors,
        verbose=verbose,
    )
    return 0


@version_option("scCountReads bins")
@bins_app.callback(invoke_without_command=True)
def bins(
    bamfiles: list[str] = typer.Option(
        ...,
        "-b",
        "--bamfiles",
        metavar=".bam",
        help="List of indexed BAM files separated by spaces.",
    ),
    barcodes: str = typer.Option(
        ...,
        "-bc",
        "--barcodes",
        metavar=".txt",
        help="A single-column file containing the cell barcode whitelist, one barcode per line.",
    ),
    outFile: str = typer.Option(
        ...,
        "-o",
        "--outFile",
        help="The file to write results to. For `scFilterStats`, `scFilterBarcodes` and `scJSD`, the output file is a .tsv file. For other tools, the output file is an updated .h5ad or .h5mu file with the result of the requested operation.",
    ),
    GTF: str = typer.Option(
        None, "--BED", metavar=".bed/.gtf", help="BED/GTF files to limit the coverage analysis to the regions in them."
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
        help="A BED or GTF file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry. Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered. Please note that you should adjust the effective genome size, if relevant.",
    ),
    chrToSkip: list[str] = typer.Option(
        None,
        "--chrToSkip",
        help="List of chromosome names to skip from the analysis. Regions on these chromosomes will be excluded. Useful for skipping mitochondrial, X chromosome, or unplaced contigs. Multiple chromosomes can be specified, e.g. ``--chrToSkip chrM chrX``.",
    ),
    binSize: int = typer.Option(
        10000, "-bs", "--binSize", help="Size of the bins, in bases, to calculate coverage. (Default: %(default)s)"
    ),
    distanceBetweenBins: int = typer.Option(
        0,
        "--distanceBetweenBins",
        help="The gap length, in bases, between bins for calculating coverage. Larger values can be used to sample the genome for input files with high coverage",
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
    minMappingQuality: int = typer.Option(
        None,
        "--minMappingQuality",
        help="If set, only reads that have a mapping quality score of at least this are considered.",
    ),
    samFlagInclude: int = typer.Option(
        None,
        "--samFlagInclude",
        help="Include reads based on SAM flag. For example, to get only reads that are the first mate, use a flag of 64. This is useful to count properly paired reads only once, as otherwise the second mate will be also considered for the coverage. This argument can be used more than once in a command. (Default: %(default)s)",
    ),
    samFlagExclude: int = typer.Option(
        None,
        "--samFlagExclude",
        help="Exclude reads based on the SAM flag. For example, to get only reads that map to the forward strand, use ``--samFlagExclude 16``, where 16 is the SAM flag for reads that map to the reverse strand. This argument can be used more than once in a command. (Default: %(default)s)",
    ),
    minFragmentLength: int = typer.Option(
        0,
        "--minFragmentLength",
        help="The minimum fragment length needed for read/pair inclusion. This option is for useful in ATACseq experiments, for filtering mono- or di-nucleosome fragments. (Default: %(default)s)",
    ),
    maxFragmentLength: int = typer.Option(
        0,
        "--maxFragmentLength",
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
        flag_value=True,
        help="This parameter allows the extension of reads to fragment size. If set, each read is extended, without exception.\n\n*NOTE*: This feature is generally NOT recommended for spliced-read data, such as RNA-seq, as it would extend reads over skipped regions.\n\n*Single-end*: Requires a user specified value for the final fragment length. Reads that already exceed this fragment length will not be extended.\n\n*Paired-end*: Reads with mates are always extended to match the fragment size defined by the two read mates. Unmated reads, mate reads that map too far apart (>4x fragment length) or even map to different chromosomes are treated like single-end reads. The input of a fragment length value is optional. If no value is specified, it is estimated from the data (mean of the fragment size of all mate reads).",
    ),
    centerReads: bool = typer.Option(
        False,
        "--centerReads",
        help="By adding this option, reads are centered with respect to the fragment length. For paired-end data, the read is centered at the fragment length defined by the two ends of the fragment. For single-end data, the given fragment length is used. This option is useful to get sharper signal around enriched regions.",
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
        help="Check whether a given motif is present in the read and the corresponding reference genome. This option checks for the motif at the 5-end of the read and at the 5-overhang in the genome, which is useful in identifying reads properly cut by a restriction-enzyme or MNAse. For example, if you want to search for an \"A\" at the 5'-end of the read and \"TA\" at 5'-overhang, use ``-m 'A,TA'``. Reads not containing the given motif are filtered out. ",
    ),
    genome2bit: str = typer.Option(
        None,
        "-g",
        "--genome2bit",
        help="If ``--motifFilter`` is provided, please also provide the genome sequence in 2bit format.",
    ),
    GCcontentFilter: str = typer.Option(
        None,
        "-gc",
        "--GCcontentFilter",
        help="Check whether the GC content of the read falls within the provided range Input format must be '<low>,<high>' , where <low> is the lower bound and <high> is the upper bound in the fraction of GC (eg. '0.1,0.9' ). If the GC content of the reads fall outside the range, they are filtered out. ",
    ),
    minAlignedFraction: float = typer.Option(
        None,
        "--minAlignedFraction",
        help="Minimum fraction of the reads which should be aligned to be counted. This includes mismatches tolerated by the aligners, but excludes InDels/Clippings. (Default: %(default)s)",
    ),
    valueTag: str = typer.Option(
        None,
        "--valueTag",
        help='Instead of counting each read/fragment as "1", add the values from a given BAM tag to the count matrix. For example, this can be used to count the number of methylated CpG per fragment.',
    ),
    genomeChunkSize: int = typer.Option(
        None,
        "--genomeChunkSize",
        help="Manually specify the size of the genome provided to each processor. The default value of None specifies that this is determined by read density of the BAM file.",
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
    return _count_reads(
        bamfiles,
        barcodes,
        outFile,
        GTF,
        labels,
        smartLabels,
        region,
        blackListFileName,
        chrToSkip,
        binSize,
        distanceBetweenBins,
        cellTag,
        groupTag,
        minMappingQuality,
        samFlagInclude,
        samFlagExclude,
        minFragmentLength,
        maxFragmentLength,
        filterRNAstrand,
        extendReads,
        centerReads,
        duplicateFilter,
        motifFilter,
        genome2bit,
        GCcontentFilter,
        minAlignedFraction,
        valueTag,
        genomeChunkSize,
        numberOfProcessors,
        verbose,
    )


@version_option("scCountReads features")
@features_app.callback(invoke_without_command=True)
def features(
    bamfiles: list[str] = typer.Option(
        ..., "-b", "--bamfiles", metavar=".bam", help="List of indexed BAM files separated by spaces."
    ),
    barcodes: str = typer.Option(
        ...,
        "-bc",
        "--barcodes",
        metavar=".txt",
        help="A single-column file containing the cell barcode whitelist, one barcode per line.",
    ),
    outFile: str = typer.Option(
        ...,
        "-o",
        "--outFile",
        help="The file to write results to. For `scFilterStats`, `scFilterBarcodes` and `scJSD`, the output file is a .tsv file. For other tools, the output file is an updated .h5ad or .h5mu file with the result of the requested operation.",
    ),
    GTF: str = typer.Option(
        ..., "--BED", metavar=".bed/.gtf", help="BED/GTF files to limit the coverage analysis to the regions in them."
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
        help="A BED or GTF file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry. Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered. Please note that you should adjust the effective genome size, if relevant.",
    ),
    chrToSkip: list[str] = typer.Option(
        None,
        "--chrToSkip",
        help="List of chromosome names to skip from the analysis. Regions on these chromosomes will be excluded. Useful for skipping mitochondrial, X chromosome, or unplaced contigs. Multiple chromosomes can be specified, e.g. ``--chrToSkip chrM chrX``.",
    ),
    binSize: int = typer.Option(
        10000, "-bs", "--binSize", help="Size of the bins, in bases, to calculate coverage. (Default: %(default)s)"
    ),
    distanceBetweenBins: int = typer.Option(
        0,
        "--distanceBetweenBins",
        help="The gap length, in bases, between bins for calculating coverage. Larger values can be used to sample the genome for input files with high coverage",
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
    minMappingQuality: int = typer.Option(
        None,
        "--minMappingQuality",
        help="If set, only reads that have a mapping quality score of at least this are considered.",
    ),
    samFlagInclude: int = typer.Option(
        None,
        "--samFlagInclude",
        help="Include reads based on SAM flag. For example, to get only reads that are the first mate, use a flag of 64. This is useful to count properly paired reads only once, as otherwise the second mate will be also considered for the coverage. This argument can be used more than once in a command. (Default: %(default)s)",
    ),
    samFlagExclude: int = typer.Option(
        None,
        "--samFlagExclude",
        help="Exclude reads based on the SAM flag. For example, to get only reads that map to the forward strand, use ``--samFlagExclude 16``, where 16 is the SAM flag for reads that map to the reverse strand. This argument can be used more than once in a command. (Default: %(default)s)",
    ),
    minFragmentLength: int = typer.Option(
        0,
        "--minFragmentLength",
        help="The minimum fragment length needed for read/pair inclusion. This option is for useful in ATACseq experiments, for filtering mono- or di-nucleosome fragments. (Default: %(default)s)",
    ),
    maxFragmentLength: int = typer.Option(
        0,
        "--maxFragmentLength",
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
        flag_value=True,
        help="This parameter allows the extension of reads to fragment size. If set, each read is extended, without exception.\n\n*NOTE*: This feature is generally NOT recommended for spliced-read data, such as RNA-seq, as it would extend reads over skipped regions.\n\n*Single-end*: Requires a user specified value for the final fragment length. Reads that already exceed this fragment length will not be extended.\n\n*Paired-end*: Reads with mates are always extended to match the fragment size defined by the two read mates. Unmated reads, mate reads that map too far apart (>4x fragment length) or even map to different chromosomes are treated like single-end reads. The input of a fragment length value is optional. If no value is specified, it is estimated from the data (mean of the fragment size of all mate reads).",
    ),
    centerReads: bool = typer.Option(
        False,
        "--centerReads",
        help="By adding this option, reads are centered with respect to the fragment length. For paired-end data, the read is centered at the fragment length defined by the two ends of the fragment. For single-end data, the given fragment length is used. This option is useful to get sharper signal around enriched regions.",
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
        help="Check whether a given motif is present in the read and the corresponding reference genome. This option checks for the motif at the 5-end of the read and at the 5-overhang in the genome, which is useful in identifying reads properly cut by a restriction-enzyme or MNAse. For example, if you want to search for an \"A\" at the 5'-end of the read and \"TA\" at 5'-overhang, use ``-m 'A,TA'``. Reads not containing the given motif are filtered out. ",
    ),
    genome2bit: str = typer.Option(
        None,
        "-g",
        "--genome2bit",
        help="If ``--motifFilter`` is provided, please also provide the genome sequence in 2bit format.",
    ),
    GCcontentFilter: str = typer.Option(
        None,
        "-gc",
        "--GCcontentFilter",
        help="Check whether the GC content of the read falls within the provided range Input format must be '<low>,<high>' , where <low> is the lower bound and <high> is the upper bound in the fraction of GC (eg. '0.1,0.9' ). If the GC content of the reads fall outside the range, they are filtered out. ",
    ),
    minAlignedFraction: float = typer.Option(
        None,
        "--minAlignedFraction",
        help="Minimum fraction of the reads which should be aligned to be counted. This includes mismatches tolerated by the aligners, but excludes InDels/Clippings. (Default: %(default)s)",
    ),
    valueTag: str = typer.Option(
        None,
        "--valueTag",
        help='Instead of counting each read/fragment as "1", add the values from a given BAM tag to the count matrix. For example, this can be used to count the number of methylated CpG per fragment.',
    ),
    genomeChunkSize: int = typer.Option(
        None,
        "--genomeChunkSize",
        help="Manually specify the size of the genome provided to each processor. The default value of None specifies that this is determined by read density of the BAM file.",
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
    return _count_reads(
        bamfiles,
        barcodes,
        outFile,
        GTF,
        labels,
        smartLabels,
        region,
        blackListFileName,
        chrToSkip,
        binSize,
        distanceBetweenBins,
        cellTag,
        groupTag,
        minMappingQuality,
        samFlagInclude,
        samFlagExclude,
        minFragmentLength,
        maxFragmentLength,
        filterRNAstrand,
        extendReads,
        centerReads,
        duplicateFilter,
        motifFilter,
        genome2bit,
        GCcontentFilter,
        minAlignedFraction,
        valueTag,
        genomeChunkSize,
        numberOfProcessors,
        verbose,
    )


if __name__ == "__main__":
    app()
