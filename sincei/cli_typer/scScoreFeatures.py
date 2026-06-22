from __future__ import annotations

import typer

from sincei import _sincei as internal

from ._common_args import (
    GTF_GFF_OPTS,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    OverlapPolicy,
    log_parameters,
    override,
    preprocess_args,
)

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
    "--overlapPolicy",
    metavar="POLICY",
    rich_help_panel="Scoring options",
    help="Policy for handling regions in the .h5ad input that only partially overlap regions in --features.\n\n"
    "One of: [bold yellow]partial[/bold yellow], [bold yellow]all[/bold yellow], [bold yellow]none[/bold yellow].",
)
CENTER_SCORES = typer.Option(
    False,
    "-cs",
    "--centerScores",
    rich_help_panel="Scoring options",
    help="If set, center and scale the scores to unit variance and zero mean.",
)
FEATURES_HELP = "Path to the BED or GTF file containing the features to use for aggregation/scoring."

# GTF/GFF options shared by both subcommands (only affect GTF/GFF inputs;
# ignored for BED).  Defaults are gene-level so GTF/GFF features are keyed by
# gene rather than by transcript.
FEATURE_TYPE = override(
    GTF_GFF_OPTS["transcript_id"],
    default="gene",
    help="The column-3 feature type keyed as a region. Defaults to ``gene`` so GTF/GFF inputs are "
    "scored per gene; set to ``transcript`` (or another type) to key on different features.",
)
NAME_ATTR = override(
    GTF_GFF_OPTS["transcript_id_tag"],
    default="gene_id",
    help="The column-9 attribute key whose value names each region (and, in ``--metagene`` mode, "
    "groups exons). Defaults to ``gene_id``.",
)
EXON_ID = GTF_GFF_OPTS["exon_id"]
METAGENE = override(
    GTF_GFF_OPTS["metagene"],
    help="Merge the exons of each gene and score over the merged exons rather than the full genomic "
    "interval. Exon rows are those whose column-3 type matches ``--exonID``, grouped by the "
    "``--transcriptIDtag`` attribute.",
)


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
        "--maxRegion",
        rich_help_panel="Activity options",
        help="Maximum region size (in kb) upstream and downstream of the genes to consider for activity calculation.",
    ),
    gene_body: bool = typer.Option(
        False,
        "--geneBody",
        rich_help_panel="Activity options",
        help="If set, the entire gene body is weighted as 1 (like the TSS) and decay starts beyond the gene body. "
        "By default the weight decay starts from the TSS.",
    ),
    normalize_gene_lengths: bool = typer.Option(
        False,
        "--normalizeGeneLengths",
        rich_help_panel="Activity options",
        help="If set, gene scores are normalized with respect to gene length in the input GTF/BED file.",
    ),
    exclude_in_range: str = typer.Option(
        None,
        "--excludeInRange",
        rich_help_panel="Activity options",
        help="Exclude regions that overlap other features from contributing to the activity score of the input "
        "genes. Options are: 'TSS' (exclude features overlapping the TSS of other genes) or 'genes' (exclude features "
        "overlapping the bodies of other genes).",
    ),
    feature_type: str = FEATURE_TYPE,
    exon_id: str = EXON_ID,
    name_attr: str = NAME_ATTR,
    metagene: bool = METAGENE,
    number_of_processors: str = OTHER_OPTS["number_of_processors"],
    verbose: bool = OTHER_OPTS["verbose"],
    help: bool = OTHER_OPTS["help"],
) -> int:
    if verbose:
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
            feature_type=feature_type,
            exon_id=exon_id,
            name_attr=name_attr,
            metagene=metagene,
            number_of_processors=number_of_processors,
        )
    # In metagene mode the kept rows are exons (--exonID); otherwise the region
    # type itself (--transcriptID, default 'gene'). The backend groups exons by
    # name_attr when metagene is set.
    region_type = exon_id if metagene else feature_type
    result_path = internal.score_features(
        input,
        features,
        out_file,
        mode="activities",
        overlap_policy=overlap_policy.value,
        center_scores=center_scores,
        decay=decay,
        max_region=max_region,
        gene_body=gene_body,
        normalize_gene_lengths=normalize_gene_lengths,
        exclude_in_range=exclude_in_range,
        feature_type=region_type,
        name_attr=name_attr,
        metagene=metagene,
        num_threads=number_of_processors,
    )
    typer.echo(result_path)
    return 0


@aggregate_app.callback(invoke_without_command=True)
def aggregate(
    input: str = INPUT_OUTPUT_OPTS["h5ad_file"],
    out_file: str = INPUT_OUTPUT_OPTS["out_file"],
    features: str = typer.Option(
        ..., "--features", metavar=".bed/.gtf/.gff", rich_help_panel="Scoring options", help=FEATURES_HELP
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
    feature_type: str = FEATURE_TYPE,
    exon_id: str = EXON_ID,
    name_attr: str = NAME_ATTR,
    metagene: bool = METAGENE,
    number_of_processors: str = OTHER_OPTS["number_of_processors"],
    verbose: bool = OTHER_OPTS["verbose"],
    help: bool = OTHER_OPTS["help"],
) -> int:
    if verbose:
        log_parameters(
            input=input,
            out_file=out_file,
            features=features,
            overlap_policy=overlap_policy,
            center_scores=center_scores,
            penalty=penalty,
            feature_type=feature_type,
            exon_id=exon_id,
            name_attr=name_attr,
            metagene=metagene,
            number_of_processors=number_of_processors,
        )
    # In metagene mode the kept rows are exons (--exonID); otherwise the region
    # type itself (--transcriptID, default 'gene'). The backend groups exons by
    # name_attr when metagene is set.
    region_type = exon_id if metagene else feature_type
    result_path = internal.score_features(
        input,
        features,
        out_file,
        mode="aggregate",
        overlap_policy=overlap_policy.value,
        center_scores=center_scores,
        penalty=penalty,
        feature_type=region_type,
        name_attr=name_attr,
        metagene=metagene,
        num_threads=number_of_processors,
    )
    typer.echo(result_path)
    return 0


def cli() -> None:
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
