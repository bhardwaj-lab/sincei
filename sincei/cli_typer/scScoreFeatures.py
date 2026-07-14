from __future__ import annotations

import typer

from sincei import _sincei as internal

from ._common_args import (
    GTF_GFF_OPTS,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    OverlapPolicy,
    log_parameters,
    preprocess_args,
)

DESCRIPTION = (
    "Aggregate a binned chromatin count matrix into per-feature scores.\n\n"
    "``scScoreFeatures`` sums the counts of the bins overlapping each feature in ``--features``, "
    "producing a cells x features matrix. The features can be genes (from a GTF/GFF) or Variable "
    "Chromatin Regions (from a BED file produced by scFindVCRs)."
)


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    context_settings={"help_option_names": []},
)


@app.command(help=DESCRIPTION)
def main(
    input: str = INPUT_OUTPUT_OPTS["h5ad_file"],
    out_file: str = INPUT_OUTPUT_OPTS["out_file"],
    features: str = typer.Option(
        ...,
        "--features",
        metavar=".bed/.gtf/.gff",
        rich_help_panel="Scoring options",
        help="Path to the BED, GTF or GFF file containing the features to score.",
    ),
    overlap_policy: OverlapPolicy = typer.Option(
        OverlapPolicy.partial,
        "-op",
        "--overlapPolicy",
        metavar="POLICY",
        rich_help_panel="Scoring options",
        help="How to treat a bin in the .h5ad input that only partially overlaps a region in --features.\n\n"
        "[bold yellow]partial[/bold yellow]: count the fraction of the bin lying inside the region; "
        "[bold yellow]all[/bold yellow]: count the whole bin; "
        "[bold yellow]none[/bold yellow]: ignore the bin unless it lies wholly inside the region.",
    ),
    penalty: float = typer.Option(
        None,
        "-pen",
        "--penalty",
        rich_help_panel="Scoring options",
        help="Penalty value to determine which VCRs to score. Used only when the input is a BED file "
        "created with ``scFindVCRs`` with a range of penalties (stored in the 5th column).",
    ),
    # GTF/GFF options; only affect GTF/GFF inputs, ignored for BED.
    feature_type: list[str] | None = GTF_GFF_OPTS["transcript_id"],
    exon_id: list[str] | None = GTF_GFF_OPTS["exon_id"],
    name_attr: str | None = GTF_GFF_OPTS["transcript_id_tag"],
    metagene: bool = GTF_GFF_OPTS["metagene"],
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
            penalty=penalty,
            feature_type=feature_type,
            exon_id=exon_id,
            name_attr=name_attr,
            metagene=metagene,
            number_of_processors=number_of_processors,
        )
    result_path = internal.score_features(
        input,
        features,
        out_file,
        overlap_policy=overlap_policy.value,
        penalty=penalty,
        feature_type=feature_type,
        exon_type=exon_id,
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
