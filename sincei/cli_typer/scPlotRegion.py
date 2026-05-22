from __future__ import annotations

from typing import Annotated

import typer

from .shared import log_parameters, normalize_processors, normalize_region, version_option

DESCRIPTION = """
Plot pseudo-bulk and per cell coverage for a genomic region.
"""


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help=DESCRIPTION,
    context_settings={"help_option_names": ["-h", "--help"]},
)


@version_option("scPlotRegion")
@app.callback(invoke_without_command=True)
def main(
    input: Annotated[
        str, typer.Option("-i", "--input", metavar=".h5ad", help="sincei-generated input file in .h5ad format.")
    ],
    outFile: Annotated[str, typer.Option("-o", "--outFile", help="The file to write results to.")],
    region: Annotated[
        str,
        typer.Option(
            "-r",
            "--region",
            callback=normalize_region,
            help="Region of the genome to limit the operation to - this is useful when testing parameters to reduce the computing time. The format is chr:start:end, for example ``--region chr10`` or ``--region chr10:456700:891000``.",
        ),
    ],
    mode: Annotated[str, typer.Option("-m", "--mode", help="Aggregation mode for the top subplot")],
    summarize: bool = typer.Option(
        False, "--summarize", help="Sum contiguous columns to reduce matrix width for heatmap plotting (default: False)"
    ),
    signal_min: float = typer.Option(
        None, "--signal-min", help="Minimum value for pseudobulk track plot (default: data min)"
    ),
    signal_max: float = typer.Option(
        None, "--signal-max", help="Maximum value for pseudobulk track plot (default: data max)"
    ),
    map_min: float = typer.Option(None, "--map-min", help="Minimum value for cell heatmap (default: data min)"),
    map_max: float = typer.Option(None, "--map-max", help="Maximum value for cell heatmap (default: data max)"),
    color: str = typer.Option("red", "--color", help="Color for top line plot (default: red)"),
    colormap: str = typer.Option("Reds", "--colormap", help="Colormap for heatmap (default: Redss)"),
    figsize: tuple[float, float] = typer.Option(
        (14, 8), "--figsize", help="Figure size in inches (width height, default: 14 8)"
    ),
    dpi: int = typer.Option(300, "--dpi", help="DPI for output PNG (default: 300)"),
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
        region=region,
        mode=mode,
        summarize=summarize,
        signal_min=signal_min,
        signal_max=signal_max,
        map_min=map_min,
        map_max=map_max,
        color=color,
        colormap=colormap,
        figsize=figsize,
        dpi=dpi,
        numberOfProcessors=numberOfProcessors,
        verbose=verbose,
    )
    return 0


if __name__ == "__main__":
    app()
