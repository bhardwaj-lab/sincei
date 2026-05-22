from __future__ import annotations

import typer

from .shared import log_parameters, normalize_processors, version_option

DESCRIPTION = """
Calculate gene/region activity scores from binned chromatin count matrix.
"""


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help=DESCRIPTION,
    context_settings={"help_option_names": ["-h", "--help"]},
)

activities_app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Calculate gene activity scores.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
aggregate_app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Aggregate binned chromatin data into VCRs.",
    context_settings={"help_option_names": ["-h", "--help"]},
)

app.add_typer(activities_app, name="activities")
app.add_typer(aggregate_app, name="aggregate")


@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    input: str = typer.Option(
        None, "-i", "--input", metavar=".h5ad", help="sincei-generated input file in .h5ad format."
    ),
    outFile: str = typer.Option(None, "-o", "--outFile", help="The file to write results to."),
    GTF: str = typer.Option(None, "--features", metavar=".bed/.gtf"),
    overlapPolicy: str = typer.Option("partial", "--overlapPolicy", "-op"),
    centerScores: bool = typer.Option(False, "--centerScores", "-cs"),
    numberOfProcessors: str = typer.Option("max", "-p", "--numberOfProcessors", callback=normalize_processors),
    verbose: bool = typer.Option(False, "-v", "--verbose"),
) -> int:
    if ctx.invoked_subcommand is None:
        if input is None and outFile is None and GTF is None:
            typer.echo(ctx.get_help())
            raise typer.Exit()
        log_parameters(
            input=input,
            outFile=outFile,
            GTF=GTF,
            overlapPolicy=overlapPolicy,
            centerScores=centerScores,
            numberOfProcessors=numberOfProcessors,
            verbose=verbose,
        )
    return 0


@version_option("scScoreFeatures activities")
@activities_app.callback(invoke_without_command=True)
def activities(
    input: str = typer.Option(
        ..., "-i", "--input", metavar=".h5ad", help="sincei-generated input file in .h5ad format."
    ),
    outFile: str = typer.Option(..., "-o", "--outFile", help="The file to write results to."),
    GTF: str = typer.Option(
        ...,
        "--features",
        metavar=".bed/.gtf",
        help="Path to the BED or GTF file containing the features to use for aggregation/scoring.",
    ),
    overlapPolicy: str = typer.Option(
        "partial",
        "--overlapPolicy",
        "-op",
        help="Policy for handling regions present in .h5ad input file that only partially overlap regions present in --features. Options are: partial, all, none. Default: %(default)s.",
    ),
    centerScores: bool = typer.Option(
        False,
        "--centerScores",
        "-cs",
        help="If flag is set, center and scale the scores to unit variance and zero mean. Default: %(default)s.",
    ),
    decay: float = typer.Option(
        0.75,
        "-d",
        "--decay",
        help="Decay parameter for calculating distance weights. Higher values lead to faster decay with distance. Weights are calculated as ``exp(-decay * distance_in_kb / 10)``. Only used with ``--mode activities``. Default: %(default)s.",
    ),
    maxRegion: int = typer.Option(
        100,
        "-mr",
        "--maxRegion",
        help="Maximum region size (in kb) upstream and downstream of the genes to consider for activity calculation. Default: %(default)s.",
    ),
    geneBody: bool = typer.Option(
        False,
        "--geneBody",
        help="Flag to indicate whether the entire gene body is weighted as 1 (like the TSS). If provided, decay starts beyond gene body. By default, the weight decay starts from TSS. Default: %(default)s.",
    ),
    normalizeGeneLengths: bool = typer.Option(
        False,
        "--normalizeGeneLengths",
        help="Flag to indicate whether to apply length normalization to the input genes. If provided, gene scores are normalized w.r.t. gene length in the input GTF/BED file. Default: %(default)s.",
    ),
    excludeInRange: str = typer.Option(
        None,
        "--excludeInRange",
        help="Exclude regions that overlap other features from contributing to activity score of the input genes. This could help avoid spurious correlations between the target genes and the neighboring genes (in particular for promoter-enriched signals, such as H3K4me3). Options are: 'TSS': exclude features overlapping the TSS of other genes. 'genes': exclude features overlapping the bodies of other genes. Default: %(default)s.",
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
        input=input,
        outFile=outFile,
        GTF=GTF,
        overlapPolicy=overlapPolicy,
        centerScores=centerScores,
        decay=decay,
        maxRegion=maxRegion,
        geneBody=geneBody,
        normalizeGeneLengths=normalizeGeneLengths,
        excludeInRange=excludeInRange,
        numberOfProcessors=numberOfProcessors,
        verbose=verbose,
    )
    return 0


@version_option("scScoreFeatures aggregate")
@aggregate_app.callback(invoke_without_command=True)
def aggregate(
    input: str = typer.Option(
        ..., "-i", "--input", metavar=".h5ad", help="sincei-generated input file in .h5ad format."
    ),
    outFile: str = typer.Option(..., "-o", "--outFile", help="The file to write results to."),
    GTF: str = typer.Option(
        ...,
        "--features",
        metavar=".bed/.gtf",
        help="Path to the BED or GTF file containing the features to use for aggregation/scoring.",
    ),
    overlapPolicy: str = typer.Option(
        "partial",
        "--overlapPolicy",
        "-op",
        help="Policy for handling regions present in .h5ad input file that only partially overlap regions present in --features. Options are: partial, all, none. Default: %(default)s.",
    ),
    centerScores: bool = typer.Option(
        False,
        "--centerScores",
        "-cs",
        help="If flag is set, center and scale the scores to unit variance and zero mean. Default: %(default)s.",
    ),
    penalty: float = typer.Option(
        None,
        "-pen",
        "--penalty",
        help="Penalty value to determine which VCRs to use for aggregation. Used only when the input is a BED file created with ``scFindVCRs`` with a range of penalties (stored in the 5th column). Default: %(default)s.",
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
        input=input,
        outFile=outFile,
        GTF=GTF,
        overlapPolicy=overlapPolicy,
        centerScores=centerScores,
        penalty=penalty,
        numberOfProcessors=numberOfProcessors,
        verbose=verbose,
    )
    return 0


if __name__ == "__main__":
    app()
