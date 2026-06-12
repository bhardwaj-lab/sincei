from __future__ import annotations

import typer

from ._common_args import INPUT_OUTPUT_OPTS, OTHER_OPTS, OverlapPolicy, log_parameters, override, preprocess_args

DESCRIPTION = (
    "Calculate gene/region activity scores from a binned chromatin count matrix.\n\n"
    "``scScoreFeatures`` computes gene activity scores from chromatin data or aggregates"
    "binned chromatin data into Variable Chromatin Regions (use output from scFindVCRs)."
)


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)

activities_app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help="Calculate gene activity scores.",
    context_settings={"help_option_names": []},
)
aggregate_app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help="Aggregate binned chromatin data into VCRs.",
    context_settings={"help_option_names": []},
)

app.add_typer(activities_app, name="activities")
app.add_typer(aggregate_app, name="aggregate")

# Options shared by both subcommands.
OVERLAP_POLICY = typer.Option(
    OverlapPolicy.partial,
    "-op",
    "--overlap-policy",
    metavar="POLICY",
    rich_help_panel="Scoring options",
    help="Policy for handling regions in the .h5ad input that only partially overlap regions in --features.\n\n"
    "One of: [bold yellow]partial[/bold yellow], [bold yellow]all[/bold yellow], [bold yellow]none[/bold yellow].",
)
CENTER_SCORES = typer.Option(
    False,
    "-cs",
    "--center-scores",
    rich_help_panel="Scoring options",
    help="If set, center and scale the scores to unit variance and zero mean.",
)
FEATURES_HELP = "Path to the BED or GTF file containing the features to use for aggregation/scoring."


@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    input: str = override(INPUT_OUTPUT_OPTS["h5ad_file"], default=None),
    out_file: str = override(INPUT_OUTPUT_OPTS["out_file"], default=None),
    features: str = typer.Option(
        None, "--features", metavar=".bed/.gtf", rich_help_panel="Scoring options", help=FEATURES_HELP
    ),
    overlap_policy: OverlapPolicy = OVERLAP_POLICY,
    center_scores: bool = CENTER_SCORES,
    number_of_processors: str = OTHER_OPTS["number_of_processors"],
    verbose: bool = OTHER_OPTS["verbose"],
    help: bool = OTHER_OPTS["help"],
) -> int:
    if ctx.invoked_subcommand is None:
        if input is None and out_file is None and features is None:
            typer.echo(ctx.get_help())
            raise typer.Exit()
        log_parameters(
            input=input,
            out_file=out_file,
            features=features,
            overlap_policy=overlap_policy,
            center_scores=center_scores,
            number_of_processors=number_of_processors,
            verbose=verbose,
        )
    return 0


@activities_app.callback(invoke_without_command=True)
def activities(
    input: str = INPUT_OUTPUT_OPTS["h5ad_file"],
    out_file: str = INPUT_OUTPUT_OPTS["out_file"],
    features: str = typer.Option(
        ..., "--features", metavar=".bed/.gtf", rich_help_panel="Scoring options", help=FEATURES_HELP
    ),
    overlap_policy: OverlapPolicy = OVERLAP_POLICY,
    center_scores: bool = CENTER_SCORES,
    decay: float = typer.Option(
        0.75,
        "-d",
        "--decay",
        rich_help_panel="Activity options",
        help="Decay parameter for calculating distance weights. Higher values lead to faster decay with distance. "
        "Weights are calculated as ``exp(-decay * distance_in_kb / 10)``.",
    ),
    max_region: int = typer.Option(
        100,
        "-mr",
        "--max-region",
        rich_help_panel="Activity options",
        help="Maximum region size (in kb) upstream and downstream of the genes to consider for activity calculation.",
    ),
    gene_body: bool = typer.Option(
        False,
        "--gene-body",
        rich_help_panel="Activity options",
        help="If set, the entire gene body is weighted as 1 (like the TSS) and decay starts beyond the gene body. "
        "By default the weight decay starts from the TSS.",
    ),
    normalize_gene_lengths: bool = typer.Option(
        False,
        "--normalize-gene-lengths",
        rich_help_panel="Activity options",
        help="If set, gene scores are normalized with respect to gene length in the input GTF/BED file.",
    ),
    exclude_in_range: str = typer.Option(
        None,
        "--exclude-in-range",
        rich_help_panel="Activity options",
        help="Exclude regions that overlap other features from contributing to the activity score of the input "
        "genes. Options are: 'TSS' (exclude features overlapping the TSS of other genes) or 'genes' (exclude features "
        "overlapping the bodies of other genes).",
    ),
    number_of_processors: str = OTHER_OPTS["number_of_processors"],
    verbose: bool = OTHER_OPTS["verbose"],
    help: bool = OTHER_OPTS["help"],
) -> int:
    log_parameters(
        input=input,
        out_file=out_file,
        features=features,
        overlap_policy=overlap_policy,
        center_scores=center_scores,
        decay=decay,
        max_region=max_region,
        gene_body=gene_body,
        normalize_gene_lengths=normalize_gene_lengths,
        exclude_in_range=exclude_in_range,
        number_of_processors=number_of_processors,
        verbose=verbose,
    )
    return 0


@aggregate_app.callback(invoke_without_command=True)
def aggregate(
    input: str = INPUT_OUTPUT_OPTS["h5ad_file"],
    out_file: str = INPUT_OUTPUT_OPTS["out_file"],
    features: str = typer.Option(
        ..., "--features", metavar=".bed/.gtf", rich_help_panel="Scoring options", help=FEATURES_HELP
    ),
    overlap_policy: OverlapPolicy = OVERLAP_POLICY,
    center_scores: bool = CENTER_SCORES,
    penalty: float = typer.Option(
        None,
        "-pen",
        "--penalty",
        rich_help_panel="Aggregate options",
        help="Penalty value to determine which VCRs to use for aggregation. Used only when the input is a BED file "
        "created with ``scFindVCRs`` with a range of penalties (stored in the 5th column).",
    ),
    number_of_processors: str = OTHER_OPTS["number_of_processors"],
    verbose: bool = OTHER_OPTS["verbose"],
    help: bool = OTHER_OPTS["help"],
) -> int:
    log_parameters(
        input=input,
        out_file=out_file,
        features=features,
        overlap_policy=overlap_policy,
        center_scores=center_scores,
        penalty=penalty,
        number_of_processors=number_of_processors,
        verbose=verbose,
    )
    return 0


def cli() -> None:
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
