from __future__ import annotations

from typing import Annotated

import typer

from ._common_args import (
    AVAILABLE_PROCESSORS,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    OverlapPolicy,
    configure_logging,
    log_parameters,
    preprocess_args,
)

DESCRIPTION = (
    "Aggregate region-level signal into gene-level scores.\n\n"
    "``scScoreFeatures`` aggregates region-level signal in a pre-processed .h5ad "
    "object into gene-level scores based on a user-provided BED/GTF file."
)


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)

_SCORING = "Common options"


@app.callback(invoke_without_command=True)
def main(
    input: Annotated[str, INPUT_OUTPUT_OPTS["h5ad_file"]],
    out_file: Annotated[str, INPUT_OUTPUT_OPTS["out_file"]],
    features: Annotated[
        str,
        typer.Option(
            "-f",
            "--features",
            metavar=".bed/.gtf",
            rich_help_panel=_SCORING,
            help=(
                "Path to the BED or GTF file containing the features to use for "
                "aggregation/scoring."
            ),
        ),
    ],
    overlap_policy: Annotated[
        OverlapPolicy,
        typer.Option(
            "-op",
            "--overlapPolicy",
            metavar="POLICY",
            rich_help_panel=_SCORING,
            help=(
                "Policy for handling regions present in the .h5ad input file that only "
                "partially overlap regions present in --features.\n\n"
                "[bold yellow]partial[/bold yellow]: count reads in anndata regions "
                "proportionally to the overlap fraction "
                "(counts_considered = feature_counts * overlap_length / "
                "region_length).\n\n"
                "[bold yellow]all[/bold yellow]: count all reads in the partially "
                "overlapping anndata regions.\n\n"
                "[bold yellow]none[/bold yellow]: only count reads in anndata regions "
                "that are fully contained within BED/GTF regions."
            ),
        ),
    ] = OverlapPolicy.partial,
    center_scores: Annotated[
        bool,
        typer.Option(
            "-cs",
            "--centerScores",
            rich_help_panel=_SCORING,
            help=(
                "If set, center and scale the scores to unit variance and zero mean."
            ),
        ),
    ] = False,
    bed_score_filter: Annotated[
        list[float] | None,
        typer.Option(
            "-bsf",
            "--bedScoreFilter",
            metavar="FLOAT",
            rich_help_panel=_SCORING,
            help=(
                "Provide a range (two values separated by space), or a threshold "
                "(upper limit) of score to determine which input features to consider "
                "for scoring. Used only when the input is a BED file containing scores "
                "(stored in the 5th column)."
            ),
        ),
    ] = None,
    max_region: Annotated[
        int,
        typer.Option(
            "-mr",
            "--maxRegion",
            metavar="INT",
            rich_help_panel=_SCORING,
            help=(
                "Maximum region size (in kb) upstream and downstream of the genes to "
                "consider for activity calculation."
            ),
        ),
    ] = 100,
    normalize_gene_lengths: Annotated[
        bool,
        typer.Option(
            "--normalizeGeneLengths",
            rich_help_panel=_SCORING,
            help=(
                "Apply length normalization to the input genes. If provided, gene "
                "scores are normalized w.r.t. gene length in the input GTF/BED file."
            ),
        ),
    ] = False,
    number_of_processors: Annotated[
        int, OTHER_OPTS["number_of_processors"]
    ] = AVAILABLE_PROCESSORS,
    verbose: Annotated[bool, OTHER_OPTS["verbose"]] = False,
    help: Annotated[bool, OTHER_OPTS["help"]] = False,
) -> int:
    log_parameters(
        input=input,
        out_file=out_file,
        features=features,
        overlap_policy=overlap_policy,
        center_scores=center_scores,
        bed_score_filter=bed_score_filter,
        max_region=max_region,
        normalize_gene_lengths=normalize_gene_lengths,
        number_of_processors=number_of_processors,
        verbose=verbose,
    )
    return 0


def cli() -> None:
    configure_logging()
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
