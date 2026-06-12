from __future__ import annotations

import typer

from ._common_args import OTHER_OPTS, preprocess_args, version_string

# Import subcommand apps to register them with the main app.
from .scBulkCoverage import app as scBulkCoverage_app
from .scClusterCells import app as scClusterCells_app
from .scCombineMods import app as scCombineMods_app
from .scCombineSamples import app as scCombineSamples_app
from .scCountQC import app as scCountQC_app
from .scCountReads import app as scCountReads_app
from .scExportSignal import app as scExportSignal_app
from .scFilterBarcodes import app as scFilterBarcodes_app
from .scFilterStats import app as scFilterStats_app
from .scFindVCRs import app as scFindVCRs_app
from .scJSD import app as scJSD_app
from .scPlotRegion import app as scPlotRegion_app
from .scReduceDims import app as scReduceDims_app
from .scScoreFeatures import app as scScoreFeatures_app

DESCRIPTION = f"""
sincei is a suite of command-line tools to explore single-cell genomics data.
Every tool name begins with the prefix `sc` followed by <tool_name>, such as:
    scBulkCoverage -b reads.bam -g groupinfo.txt -o coverage

    sincei {version_string()}
"""


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)

# Register subcommand apps.
# The order of subcommands in the help message follows the order of registration.
app.add_typer(scFilterBarcodes_app, name="scFilterBarcodes")
app.add_typer(scFilterStats_app, name="scFilterStats")
app.add_typer(scJSD_app, name="scJSD")
app.add_typer(scCountReads_app, name="scCountReads")
app.add_typer(scCountQC_app, name="scCountQC")
app.add_typer(scFindVCRs_app, name="scFindVCRs")
app.add_typer(scReduceDims_app, name="scReduceDims")
app.add_typer(scClusterCells_app, name="scClusterCells")
app.add_typer(scScoreFeatures_app, name="scScoreFeatures")
app.add_typer(scCombineSamples_app, name="scCombineSamples")
app.add_typer(scCombineMods_app, name="scCombineMods")
app.add_typer(scPlotRegion_app, name="scPlotRegion")
app.add_typer(scBulkCoverage_app, name="scBulkCoverage")
app.add_typer(scExportSignal_app, name="scExportSignal")


@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    version: bool = OTHER_OPTS["version"],
    help: bool = OTHER_OPTS["help"],
) -> int:
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
        raise typer.Exit()
    return 0


def cli() -> None:
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
