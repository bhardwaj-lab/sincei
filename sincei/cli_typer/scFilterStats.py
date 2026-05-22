from __future__ import annotations

import typer

from .shared import log_parameters, normalize_processors, version_option

DESCRIPTION = """
Produce per-cell statistics after filtering reads by user-defined criteria.
"""


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help=DESCRIPTION,
    context_settings={"help_option_names": ["-h", "--help"]},
)


@version_option("scFilterStats")
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
    binSize: int = typer.Option(100000, "-bs", "--binSize", help="Size of the bins, in bases, to calculate coverage."),
    distanceBetweenBins: int = typer.Option(1000000, "--distanceBetweenBins", help="Gap length between bins."),
    duplicateFilter: str = typer.Option(None, "--duplicateFilter", help="How to filter for duplicates?"),
    motifFilter: list[str] = typer.Option(
        None, "-m", "--motifFilter", metavar="STR", help="Motif-based filtering options."
    ),
    genome2bit: str = typer.Option(
        None, "-g", "--genome2bit", metavar=".2bit", help="Genome 2bit file used for motif filtering."
    ),
    GCcontentFilter: str = typer.Option(
        None, "-gc", "--GCcontentFilter", metavar="STR", help="GC-content filtering options."
    ),
    minAlignedFraction: float = typer.Option(
        None,
        "--minAlignedFraction",
        metavar="FLOAT",
        help="Minimum aligned fraction required for a read to be considered.",
    ),
    minMappingQuality: int = typer.Option(
        None,
        "--minMappingQuality",
        metavar="INT",
        help="If set, only reads that have a mapping quality score of at least this are considered.",
    ),
    samFlagInclude: int = typer.Option(None, "--samFlagInclude", metavar="INT", help="SAM flags that must be present."),
    samFlagExclude: int = typer.Option(None, "--samFlagExclude", metavar="INT", help="SAM flags that must be absent."),
    minFragmentLength: int = typer.Option(0, "--minFragmentLength", metavar="INT", help="Minimum fragment length."),
    maxFragmentLength: int = typer.Option(0, "--maxFragmentLength", metavar="INT", help="Maximum fragment length."),
    filterRNAstrand: str = typer.Option(
        None, "--filterRNAstrand", help="Strand filtering options for RNA-derived data."
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
        help="Region of the genome to limit the operation to - this is useful when testing parameters to reduce the computing time. The format is chr:start:end, for example ``--region chr10`` or ``--region chr10:456700:891000``.",
    ),
    blackListFileName: list[str] = typer.Option(
        None,
        "-bl",
        "--blacklist",
        "--blackListFileName",
        help="BED/GTF files to limit analysis to the provided regions.",
    ),
    chrToSkip: list[str] = typer.Option(
        None, "--chrToSkip", help="A space separated list of chromosomes to exclude. eg. chrM chrUn"
    ),
    extendReads: int = typer.Option(
        False, "-e", "--extendReads", flag_value=True, help="Extend reads before counting."
    ),
    centerReads: bool = typer.Option(False, "--centerReads", help="Center reads before counting."),
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
        labels=labels,
        smartLabels=smartLabels,
        region=region,
        blackListFileName=blackListFileName,
        chrToSkip=chrToSkip,
        extendReads=extendReads,
        centerReads=centerReads,
        numberOfProcessors=numberOfProcessors,
        verbose=verbose,
    )
    return 0


if __name__ == "__main__":
    app()
