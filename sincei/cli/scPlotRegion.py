from __future__ import annotations

import sys
from typing import Annotated

import anndata as ad
import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import typer

from sincei.plotting._plot_region import plot_region

from ._common_args import (
    CM_PER_INCH,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    PLOT_OPTS,
    PlotFileFormat,
    SummaryMode,
    configure_logging,
    log_parameters,
    override,
    preprocess_args,
)
from ._parsers import validate_anndata

DESCRIPTION = (
    "Plot pseudo-bulk and per cell coverage for a genomic region.\n\n"
    "``scPlotRegion`` plots the signal of individual cells in a genomic region as a "
    "heatmap, in the style of a track plot, together with the summary profile "
    "(pseudo-bulk signal) on top."
)


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)

_DISPLAY = "Display options"
_COLOR = "Color / Scale options"


@app.callback(invoke_without_command=True)
def main(
    input: Annotated[str, INPUT_OUTPUT_OPTS["h5ad_file"]],
    out_file: Annotated[str, INPUT_OUTPUT_OPTS["out_file"]],
    region: Annotated[str, INPUT_OUTPUT_OPTS["region"]],
    mode: Annotated[
        SummaryMode,
        typer.Option(
            "-m",
            "--mode",
            metavar="MODE",
            rich_help_panel=_DISPLAY,
            help=(
                "How to aggregate the signal of the cells for the summary profile on "
                "top.\n\n"
                "One of: [bold yellow]sum[/bold yellow], "
                "[bold yellow]mean[/bold yellow]."
            ),
        ),
    ] = SummaryMode.sum,
    signal_min: Annotated[
        float | None,
        typer.Option(
            "--signalMin",
            metavar="FLOAT",
            rich_help_panel=_COLOR,
            show_default="minimum signal in the region",
            help="Minimum value for the summary profile.",
        ),
    ] = None,
    signal_max: Annotated[
        float | None,
        typer.Option(
            "--signalMax",
            metavar="FLOAT",
            rich_help_panel=_COLOR,
            show_default="maximum signal in the region",
            help="Maximum value for the summary profile.",
        ),
    ] = None,
    map_min: Annotated[
        float | None,
        typer.Option(
            "--mapMin",
            metavar="FLOAT",
            rich_help_panel=_COLOR,
            show_default="minimum signal in the region",
            help="Minimum value for the single-cell heatmap.",
        ),
    ] = None,
    map_max: Annotated[
        float | None,
        typer.Option(
            "--mapMax",
            metavar="FLOAT",
            rich_help_panel=_COLOR,
            show_default="maximum signal in the region",
            help="Maximum value for the single-cell heatmap.",
        ),
    ] = None,
    color: Annotated[
        str,
        typer.Option(
            "--color",
            metavar="STR",
            rich_help_panel=_COLOR,
            help="Color for the summary profile.",
        ),
    ] = "red",
    colormap: Annotated[
        str,
        typer.Option(
            "--colormap",
            metavar="STR",
            rich_help_panel=_COLOR,
            help="Colormap for the heatmap. Must be a valid matplotlib colormap.",
        ),
    ] = "Reds",
    plot_width: Annotated[float, PLOT_OPTS["plot_width"]] = 36,
    plot_height: Annotated[float, PLOT_OPTS["plot_height"]] = 20,
    plot_file_format: Annotated[
        PlotFileFormat | None,
        override(
            PLOT_OPTS["plot_file_format"],
            show_default="inferred from the --outFile suffix",
            help=(
                "Image format type. If given, this option overrides the image format "
                "inferred from the suffix of --outFile.\n\n"
                "One of: [bold yellow]png[/bold yellow], "
                "[bold yellow]jpg[/bold yellow], "
                "[bold yellow]svg[/bold yellow], [bold yellow]pdf[/bold yellow]."
            ),
        ),
    ] = None,
    dpi: Annotated[int, PLOT_OPTS["dpi"]] = 300,
    verbose: Annotated[bool, OTHER_OPTS["verbose"]] = False,
    help: Annotated[bool, OTHER_OPTS["help"]] = False,
) -> int:
    if verbose:
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
            plot_width=plot_width,
            plot_height=plot_height,
            plot_file_format=plot_file_format,
            dpi=dpi,
        )

    adata = validate_anndata(ad.read_h5ad(input), input)

    try:
        figure = plot_region(
            adata,
            region,
            mode=mode.value,
            color=color,
            colormap=colormap,
            signal_min=signal_min,
            signal_max=signal_max,
            map_min=map_min,
            map_max=map_max,
            figsize=(plot_width / CM_PER_INCH, plot_height / CM_PER_INCH),
        )
    except ValueError as exc:
        sys.stderr.write(f"{exc}\n")
        raise typer.Exit(code=1) from exc

    figure.savefig(
        out_file,
        dpi=dpi,
        bbox_inches="tight",
        format=plot_file_format.value if plot_file_format else None,
    )
    plt.close(figure)
    return 0


def cli() -> None:
    configure_logging()
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
