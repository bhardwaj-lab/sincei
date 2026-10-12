from __future__ import annotations

from sincei import _sincei as internal

from . import plotting as pl
from . import tools as tl

__version__ = internal.version()

__all__ = ["pl", "tl"]
