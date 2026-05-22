from __future__ import annotations

import typer

from .shared import log_parameters, normalize_processors, version_option

DESCRIPTION = """
Compared read coverages on sampled regions using the Jensen-Shannon distance.
"""


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help=DESCRIPTION,
    context_settings={"help_option_names": ["-h", "--help"]},
)


@version_option("scJSD")
@app.callback(invoke_without_command=True)
def main(
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
    binSize: int = typer.Option(10000, "-bs", "--binSize", help="Size of the bins, in bases, to calculate coverage."),
    duplicateFilter: str = typer.Option(
        None,
        "--duplicateFilter",
        help="How to filter for duplicates? Different combinations (using start/end/umi) are possible. Read start position and read barcode are always considered. Default (None) considers all reads as passing the filter. Note that in case of paired end data, both reads in the fragment are considered (and kept). So if you wish to keep only read1, combine this option with `--samFlagInclude`. ",
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
    minAlignedFraction: float = typer.Option(
        None,
        "--minAlignedFraction",
        help="Minimum fraction of the reads which should be aligned to be counted. This includes mismatches tolerated by the aligners, but excludes InDels/Clippings. (Default: %(default)s)",
    ),
    numberOfSamples: int = typer.Option(
        100000,
        "-n",
        "--numberOfSamples",
        help="The number of bins that are sampled from the genome, for which the overlapping number of reads is computed. (Default: %(default)s)",
    ),
    skipZeros: bool = typer.Option(
        False,
        "--skipZeros",
        help="If set, then regions with zero overlapping reads for *all* given BAM files are ignored. This will result in a reduced number of read counts than that specified in --numberOfSamples",
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
    log_parameters(
        bamfiles=bamfiles,
        barcodes=barcodes,
        outFile=outFile,
        cellTag=cellTag,
        groupTag=groupTag,
        labels=labels,
        smartLabels=smartLabels,
        blackListFileName=blackListFileName,
        chrToSkip=chrToSkip,
        binSize=binSize,
        duplicateFilter=duplicateFilter,
        minMappingQuality=minMappingQuality,
        samFlagInclude=samFlagInclude,
        samFlagExclude=samFlagExclude,
        minFragmentLength=minFragmentLength,
        maxFragmentLength=maxFragmentLength,
        minAlignedFraction=minAlignedFraction,
        numberOfSamples=numberOfSamples,
        skipZeros=skipZeros,
        numberOfProcessors=numberOfProcessors,
        verbose=verbose,
    )
    return 0


if __name__ == "__main__":
    app()
