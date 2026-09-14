from __future__ import annotations

import importlib
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from types import ModuleType

    from . import (
        GLMPCA,
        ExponentialFamily,
        FeatureScorer,
        FragmentFFT,
        GetStats,
        MultiModalClustering,
        ReadCounter,
        RegionQuery,
        TopicModels,
        Utilities,
        VCRfinder,
        WriteBedGraph,
    )

__all__ = [
    "GLMPCA",
    "ExponentialFamily",
    "FeatureScorer",
    "FragmentFFT",
    "GetStats",
    "MultiModalClustering",
    "ReadCounter",
    "RegionQuery",
    "TopicModels",
    "Utilities",
    "VCRfinder",
    "WriteBedGraph",
]


def __getattr__(name: str) -> ModuleType:
    """Import the tool module ``name`` on first access (PEP 562)."""
    if name not in __all__:
        msg = f"module {__name__!r} has no attribute {name!r}"
        raise AttributeError(msg)

    return importlib.import_module(f".{name}", __name__)


def __dir__() -> list[str]:
    return sorted(__all__)
