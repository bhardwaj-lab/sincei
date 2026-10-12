#!/usr/bin/env python
"""Barcode rank ("knee") plot for ``scFilterBarcodes``.

Plots per-barcode counts (the number of non-zero bins a barcode was detected
in) against their descending rank, on a log10 y-axis.  This mirrors the
"knee plot" used to separate real cells from background barcodes.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np

if TYPE_CHECKING:
    from collections.abc import Iterable

import matplotlib as mpl

mpl.use("Agg")


def plot_barcode_rank(
    barcode_counts: Iterable[tuple[str, int]],
    out_file: str,
    *,
    min_count: int | None = None,
    figsize: tuple[float, float] = (7, 6),
    dpi: int = 200,
) -> None:
    """Write a barcode rank plot to ``out_file``.

    Parameters
    ----------
    barcode_counts
        Iterable of ``(barcode, count)`` pairs.  Order is irrelevant; counts are
        sorted here.  ``count`` is the number of non-zero bins per barcode.
    out_file
        Output image path.  The format is inferred from the extension.
    min_count
        If given (and > 0), draw a red horizontal threshold line at
        ``log10(min_count)``, the cutoff used to select barcodes.
    figsize, dpi
        Standard matplotlib figure controls.
    """
    counts = np.array([c for _, c in barcode_counts], dtype=float)
    counts = counts[counts > 0]
    if counts.size == 0:
        msg = "no barcodes with non-zero counts to plot"
        raise ValueError(msg)

    counts = np.sort(counts)[::-1]
    ranks = np.arange(1, counts.size + 1)
    log_counts = np.log10(counts)

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(
        ranks, log_counts, linestyle="none", marker="o", color="grey", markersize=1.5
    )
    ax.set_xlabel("Barcode Rank", fontsize=12)
    ax.set_ylabel("No. of nonzero bins (log10)", fontsize=12)
    ax.set_title("Ranked counts (#bins) for detected Barcodes", fontsize=13)

    if min_count and min_count > 0:
        ax.axhline(np.log10(min_count), color="r", linewidth=1)

    fig.tight_layout()
    fig.savefig(out_file, dpi=dpi, pad_inches=0.2)
    plt.close(fig)
