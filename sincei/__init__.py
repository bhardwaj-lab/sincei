from __future__ import annotations

from sys import stderr
import logging

from . import preprocessing as pp
from . import tools as tl
from . import plotting as pl
from . import export as ex

from sincei import _sincei as internal

__version__ = internal.version()

__all__ = ["pp", "tl", "pl", "ex"]

logging.basicConfig(
    stream=stderr,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
)
