from __future__ import annotations

import typer

from ._common_args import INPUT_OUTPUT_OPTS, OTHER_OPTS, log_parameters, override, preprocess_args

DESCRIPTION = "Plot pseudo-bulk and per cell coverage for a genomic region."


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)


@app.callback(invoke_without_command=True)
def main(
    input: str = INPUT_OUTPUT_OPTS["h5ad_file"],
    out_file: str = INPUT_OUTPUT_OPTS["out_file"],
    region: str = override(INPUT_OUTPUT_OPTS["region"], default=...),
    mode: str = typer.Option(
        ...,
        "-m",
        "--mode",
        rich_help_panel="Display options",
        help="Aggregation mode for the top subplot.",
    ),
    summarize: bool = typer.Option(
        False,
        "--summarize",
        rich_help_panel="Display options",
        help="Sum contiguous columns to reduce matrix width for heatmap plotting.",
    ),
    signal_min: float = typer.Option(
        None,
        "--signal-min",
        rich_help_panel="Color / Scale options",
        help="Minimum value for pseudobulk track plot (default: data min).",
    ),
    signal_max: float = typer.Option(
        None,
        "--signal-max",
        rich_help_panel="Color / Scale options",
        help="Maximum value for pseudobulk track plot (default: data max).",
    ),
    map_min: float = typer.Option(
        None,
        "--map-min",
        rich_help_panel="Color / Scale options",
        help="Minimum value for cell heatmap (default: data min).",
    ),
    map_max: float = typer.Option(
        None,
        "--map-max",
        rich_help_panel="Color / Scale options",
        help="Maximum value for cell heatmap (default: data max).",
    ),
    color: str = typer.Option(
        "red",
        "--color",
        rich_help_panel="Color / Scale options",
        help="Color for the top line plot.",
    ),
    colormap: str = typer.Option(
        "Reds",
        "--colormap",
        rich_help_panel="Color / Scale options",
        help="Colormap for the heatmap.",
    ),
    fig_width: float = typer.Option(
        14.0,
        "--fig-width",
        rich_help_panel="Figure options",
        help="Figure width in inches.",
    ),
    fig_height: float = typer.Option(
        8.0,
        "--fig-height",
        rich_help_panel="Figure options",
        help="Figure height in inches.",
    ),
    dpi: int = typer.Option(
        300,
        "--dpi",
        rich_help_panel="Figure options",
        help="DPI for the output PNG.",
    ),
    number_of_processors: str = OTHER_OPTS["number_of_processors"],
    verbose: bool = OTHER_OPTS["verbose"],
    help: bool = OTHER_OPTS["help"],
) -> int:
    log_parameters(
        input=input,
        out_file=out_file,
        region=region,
        mode=mode,
        summarize=summarize,
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
