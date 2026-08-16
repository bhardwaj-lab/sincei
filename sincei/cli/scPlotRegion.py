from __future__ import annotations

from typing import Annotated

import typer

from ._common_args import (
    AVAILABLE_PROCESSORS,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    log_parameters,
    preprocess_args,
)

DESCRIPTION = "Plot pseudo-bulk and per cell coverage for a genomic region."


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)

_DISPLAY = "Display options"
_COLOR = "Color / Scale options"
_FIGURE = "Figure options"


@app.callback(invoke_without_command=True)
def main(
    input: Annotated[str, INPUT_OUTPUT_OPTS["h5ad_file"]],
    out_file: Annotated[str, INPUT_OUTPUT_OPTS["out_file"]],
    region: Annotated[str, INPUT_OUTPUT_OPTS["region"]],
    mode: Annotated[
        str,
        typer.Option(
            "-m",
            "--mode",
            rich_help_panel=_DISPLAY,
            help="Aggregation mode for the top subplot.",
        ),
    ],
    signal_min: Annotated[
        float | None,
        typer.Option(
            "--signalMin",
            rich_help_panel=_COLOR,
            help="Minimum value for pseudobulk track plot (default: data min).",
        ),
    ] = None,
    signal_max: Annotated[
        float | None,
        typer.Option(
            "--signalMax",
            rich_help_panel=_COLOR,
            help="Maximum value for pseudobulk track plot (default: data max).",
        ),
    ] = None,
    map_min: Annotated[
        float | None,
        typer.Option(
            "--mapMin",
            rich_help_panel=_COLOR,
            help="Minimum value for cell heatmap (default: data min).",
        ),
    ] = None,
    map_max: Annotated[
        float | None,
        typer.Option(
            "--mapMax",
            rich_help_panel=_COLOR,
            help="Maximum value for cell heatmap (default: data max).",
        ),
    ] = None,
    color: Annotated[
        str,
        typer.Option(
            "--color",
            rich_help_panel=_COLOR,
            help="Color for the top line plot.",
        ),
    ] = "red",
    colormap: Annotated[
        str,
        typer.Option(
            "--colormap",
            rich_help_panel=_COLOR,
            help="Colormap for the heatmap.",
        ),
    ] = "Reds",
    fig_width: Annotated[
        float,
        typer.Option(
            "--figWidth",
            rich_help_panel=_FIGURE,
            help="Figure width in inches.",
        ),
    ] = 14.0,
    fig_height: Annotated[
        float,
        typer.Option(
            "--figHeight",
            rich_help_panel=_FIGURE,
            help="Figure height in inches.",
        ),
    ] = 8.0,
    dpi: Annotated[
        int,
        typer.Option(
            "--dpi",
            rich_help_panel=_FIGURE,
            help="DPI for the output PNG.",
        ),
    ] = 300,
    number_of_processors: Annotated[
        int, OTHER_OPTS["number_of_processors"]
    ] = AVAILABLE_PROCESSORS,
    verbose: Annotated[bool, OTHER_OPTS["verbose"]] = False,
    help: Annotated[bool, OTHER_OPTS["help"]] = False,
) -> int:
    log_parameters(
        input=input,
        out_file=out_file,
        region=region,
        mode=mode,
        signal_min=signal_min,
        signal_max=signal_max,
        map_min=map_min,
        map_max=map_max,
        color=color,
        colormap=colormap,
        fig_width=fig_width,
        fig_height=fig_height,
        dpi=dpi,
        number_of_processors=number_of_processors,
        verbose=verbose,
    )
    return 0


def cli() -> None:
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
