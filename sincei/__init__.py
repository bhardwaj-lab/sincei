from __future__ import annotations

import logging
from sys import stderr

from sincei import _sincei as internal

from . import plotting as pl

# from . import export as ex
# from . import preprocessing as pp
# from . import tools as tl

__version__ = internal.version()

__all__ = ["ex", "pl", "pp", "tl"]

logging.basicConfig(
    stream=stderr,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
)
