"""Make the shared test helpers importable.

``pyproject.toml`` sets ``--import-mode=importlib``, which does not put the test
directory on ``sys.path``, so ``import _cli_testing`` would fail during
collection.  pytest always loads this file before collecting test modules, so
adding the directory here is enough.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
