from __future__ import annotations

import typer

from .shared import log_parameters, normalize_processors, version_option

DESCRIPTION = """
Filter cell barcodes from BAM file (for droplet-based single-cell seq).
"""


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help=DESCRIPTION,
    context_settings={"help_option_names": ["-h", "--help"]},
)


@version_option("scFilterBarcodes")
@app.callback(invoke_without_command=True)
def main(
    bamfile: str = typer.Option(..., "-b", "--bamfile", help="Indexed BAM file."),
    whitelist: str = typer.Option(
        None,
        "-w",
        "--whitelist",
        metavar=".txt",
        help="A single-column file containing the whitelist of barcodes to be used.",
    ),
    outFile: str = typer.Option(
        ...,
        "-o",
        "--outFile",
        help="The file to write results to. For `scFilterStats`, `scFilterBarcodes` and `scJSD`, the output file is a .tsv file. For other tools, the output file is an updated .h5ad or .h5mu file with the result of the requested operation.",
    ),
    minHammingDist: int = typer.Option(
        0,
        "-d",
        "--minHammingDist",
        metavar="INT",
        help="Minimum hamming distance to match the barcode in whitelist. Note that increasing the hamming distance really slows down the barcode detection process.",
    ),
    minCount: int = typer.Option(
        0,
        "-mc",
        "--minCount",
        metavar="INT",
        help="Minimum no. of bins with non-zero counts, in order to report a barcode. Note that this number would range from 0 to genome size/binSize. ",
    ),
    minMappingQuality: int = typer.Option(
        None,
        "-mq",
        "--minMappingQuality",
        metavar="INT",
        help="If set, only reads that have a mapping quality score of at least this are considered.",
    ),
    rankPlot: str = typer.Option(
        None,
        "-rp",
        "--rankPlot",
        metavar="STR",
        help='The output file name to plot the ranked counts per barcode (similar to the "knee plot",but counts are be the number of non-zero bins in this case).',
    ),
    cellTag: str = typer.Option(
        "BC", "-ct", "--cellTag", help="Name of the BAM tag from which to extract barcodes. (Default: %(default)s)"
    ),
    groupTag: str = typer.Option(None, "-gt", "--groupTag", help="Group-tag used in the BAM file."),
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
    binSize: int = typer.Option(100000, "-bs", "--binSize", help="Size of the bins, in bases, to calculate coverage."),
    distanceBetweenBins: int = typer.Option(None, "--distanceBetweenBins", help="Gap length between bins."),
    duplicateFilter: str = typer.Option(None, "--duplicateFilter", help="How to filter for duplicates?"),
    motifFilter: list[str] = typer.Option(
        None, "-m", "--motifFilter", metavar="STR", help="Motif-based filtering options."
    ),
    genome2bit: str = typer.Option(None, "-g", "--genome2bit", help="Genome 2bit file used for motif filtering."),
    GCcontentFilter: str = typer.Option(
        None, "-gc", "--GCcontentFilter", metavar="STR", help="GC-content filtering options."
    ),
    minAlignedFraction: float = typer.Option(
        None, "--minAlignedFraction", help="Minimum aligned fraction required for a read to be considered."
    ),
    minFragmentLength: int = typer.Option(0, "--minFragmentLength", help="Minimum fragment length."),
    maxFragmentLength: int = typer.Option(0, "--maxFragmentLength", help="Maximum fragment length."),
    filterRNAstrand: str = typer.Option(
        None, "--filterRNAstrand", help="Strand filtering options for RNA-derived data."
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
        bamfile=bamfile,
        whitelist=whitelist,
        outFile=outFile,
        minHammingDist=minHammingDist,
        minCount=minCount,
        minMappingQuality=minMappingQuality,
        rankPlot=rankPlot,
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
