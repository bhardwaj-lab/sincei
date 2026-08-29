"""The typer command-line apps, one per sincei command.

Every console script in ``pyproject.toml`` points at ``sincei.cli.<command>``,
so this package runs on every command.
"""

from __future__ import annotations

import importlib
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import typer

    from .main import app as main_app
    from .scBulkCoverage import app as scBulkCoverage_app
    from .scClusterCells import app as scClusterCells_app
    from .scCombineCounts import app as scCombineCounts_app
    from .scCountQC import app as scCountQC_app
    from .scCountReads import app as scCountReads_app
    from .scExportSignal import app as scExportSignal_app
    from .scFilterBarcodes import app as scFilterBarcodes_app
    from .scFilterStats import app as scFilterStats_app
    from .scFindVCRs import app as scFindVCRs_app
    from .scJSD import app as scJSD_app
    from .scPlotRegion import app as scPlotRegion_app
    from .scScoreFeatures import app as scScoreFeatures_app

# Exported name -> the module whose ``app`` it is.
_APPS: dict[str, str] = {
    "main_app": "main",
    "scBulkCoverage_app": "scBulkCoverage",
    "scClusterCells_app": "scClusterCells",
    "scCombineCounts_app": "scCombineCounts",
    "scCountQC_app": "scCountQC",
    "scCountReads_app": "scCountReads",
    "scExportSignal_app": "scExportSignal",
    "scFilterBarcodes_app": "scFilterBarcodes",
    "scFilterStats_app": "scFilterStats",
    "scFindVCRs_app": "scFindVCRs",
    "scJSD_app": "scJSD",
    "scPlotRegion_app": "scPlotRegion",
    "scScoreFeatures_app": "scScoreFeatures",
}

__all__ = [
    "main_app",
    "scBulkCoverage_app",
    "scClusterCells_app",
    "scCombineCounts_app",
    "scCountQC_app",
    "scCountReads_app",
    "scExportSignal_app",
    "scFilterBarcodes_app",
    "scFilterStats_app",
    "scFindVCRs_app",
    "scJSD_app",
    "scPlotRegion_app",
    "scScoreFeatures_app",
]


def __getattr__(name: str) -> typer.Typer:
    """Import the command module behind ``name`` on first access (PEP 562)."""
    module_name = _APPS.get(name)
    if module_name is None:
        msg = f"module {__name__!r} has no attribute {name!r}"
        raise AttributeError(msg)

    app = importlib.import_module(f".{module_name}", __name__).app
    # Cache it in the module namespace, so __getattr__ runs once per name.
    globals()[name] = app
    return app


def __dir__() -> list[str]:
    return sorted(__all__)
