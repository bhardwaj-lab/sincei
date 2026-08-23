from __future__ import annotations

from typing import Annotated

import typer

from sincei import _sincei as internal

from ._common_args import (
    AVAILABLE_PROCESSORS,
    GTF_GFF_OPTS,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    OverlapPolicy,
    configure_logging,
    log_parameters,
    preprocess_args,
)

DESCRIPTION = (
    "Aggregate a binned chromatin count matrix into per-feature scores.\n\n"
    "``scScoreFeatures`` sums the counts of the bins overlapping each feature in "
    "``--features``, producing a cells x features matrix. The features can be genes "
    "(from a GTF/GFF) or Variable Chromatin Regions (from a BED file produced by "
    "scFindVCRs)."
)


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    context_settings={"help_option_names": []},
)

_SCORING = "Scoring options"


@app.command(help=DESCRIPTION)
def main(
    input: Annotated[str, INPUT_OUTPUT_OPTS["h5ad_file"]],
    out_file: Annotated[str, INPUT_OUTPUT_OPTS["out_file"]],
    features: Annotated[
        str,
        typer.Option(
            "--features",
            metavar=".bed/.gtf/.gff",
            rich_help_panel=_SCORING,
            help="Path to the BED, GTF or GFF file containing the features to score.",
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
                "How to treat a bin in the .h5ad input that only partially overlaps a "
                "region in --features.\n\n"
                "[bold yellow]partial[/bold yellow]: count the fraction of the bin "
                "lying inside the region; "
                "[bold yellow]all[/bold yellow]: count the whole bin; "
                "[bold yellow]none[/bold yellow]: ignore the bin unless it lies wholly "
                "inside the region."
            ),
        ),
    ] = OverlapPolicy.partial,
    penalty: Annotated[
        float | None,
        typer.Option(
            "-pen",
            "--penalty",
            rich_help_panel=_SCORING,
            help=(
                "Penalty value to determine which VCRs to score. Used only when the "
                "input is a BED file created with ``scFindVCRs`` with a range of "
                "penalties (stored in the 5th column)."
            ),
        ),
    ] = None,
    # GTF/GFF options; only affect GTF/GFF inputs, ignored for BED.
    feature_type: Annotated[list[str] | None, GTF_GFF_OPTS["transcript_id"]] = None,
    exon_id: Annotated[list[str] | None, GTF_GFF_OPTS["exon_id"]] = None,
    name_attr: Annotated[str | None, GTF_GFF_OPTS["transcript_id_tag"]] = None,
    metagene: Annotated[bool, GTF_GFF_OPTS["metagene"]] = False,
    number_of_processors: Annotated[
        int, OTHER_OPTS["number_of_processors"]
    ] = AVAILABLE_PROCESSORS,
    verbose: Annotated[bool, OTHER_OPTS["verbose"]] = False,
    help: Annotated[bool, OTHER_OPTS["help"]] = False,
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
    configure_logging()
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
