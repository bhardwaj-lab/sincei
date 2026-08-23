#!/usr/bin/env python
"""Tests for ``sincei.pl.plot_region``.

The function is a plotting routine, so these check the parts that carry
meaning rather than the picture: which columns are selected and in what order,
that the aggregation matches a dense computation, and that the axis limits
follow the data and the user's arguments.

Run the tests::

    pytest tests/test_plot_region.py
"""

from __future__ import annotations

import matplotlib as mpl

mpl.use("Agg")

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from sincei.plotting._plot_region import plot_region


def make_adata(counts: np.ndarray, chroms: list[str], starts: list[int]) -> ad.AnnData:
    """A count matrix whose ``var`` carries the columns plot_region needs."""
    var = pd.DataFrame(
        {
            "chrom": chroms,
            "start": starts,
            "end": [s + 100 for s in starts],
        },
        index=[f"{c}:{s}-{s + 100}" for c, s in zip(chroms, starts, strict=True)],
    )
    return ad.AnnData(X=sp.csr_matrix(counts.astype("float32")), var=var)


@pytest.fixture
def adata() -> ad.AnnData:
    # Deliberately unsorted starts, so a test can tell selection order from
    # storage order. chr2 exists only to be excluded.
    counts = np.array([
        [1.0, 2.0, 3.0, 40.0],
        [4.0, 0.0, 0.0, 50.0],
        [0.0, 9.0, 0.0, 60.0],
    ])
    return make_adata(counts, ["chr1", "chr1", "chr1", "chr2"], [200, 0, 100, 0])


def axes_of(fig: mpl.figure.Figure) -> tuple:
    """The track axes and the heatmap image of the current figure."""
    return fig.axes[0], fig.axes[1].images[0]


def test_only_the_features_in_the_region_are_plotted(adata: ad.AnnData) -> None:
    plot_region(adata, "chr1")
    track, image = axes_of(plt.gcf())
    # Three chr1 bins, and the chr2 bin's large counts are absent.
    assert image.get_array().shape == (3, 3)
    assert track.get_xlim() == (0.0, 2.0)
    plt.close("all")


def test_features_are_ordered_by_coordinate_not_by_storage(adata: ad.AnnData) -> None:
    plot_region(adata, "chr1")
    track = plt.gcf().axes[0]
    # Stored as start 200, 0, 100. Sorted, the column sums are 2+9=11, 3, 1+4=5.
    line = track.lines[0]
    assert list(line.get_ydata()) == [11.0, 3.0, 5.0]
    plt.close("all")


def test_the_aggregate_matches_a_dense_computation(adata: ad.AnnData) -> None:
    for mode, reduce in (("sum", np.sum), ("mean", np.mean)):
        plot_region(adata, "chr1:0-300", mode=mode)
        line = plt.gcf().axes[0].lines[0]
        dense = np.asarray(
            adata[:, ["chr1:0-100", "chr1:100-200", "chr1:200-300"]].X.todense()
        )
        assert np.allclose(line.get_ydata(), reduce(dense, axis=0))
        plt.close("all")


def test_the_rows_are_ordered_by_their_totals(adata: ad.AnnData) -> None:
    plot_region(adata, "chr1")
    _, image = axes_of(plt.gcf())
    # Row totals over chr1 are 6, 4 and 9, so the middle row sinks to the
    # bottom. The leftmost column then reads 9, 2, 0 from the top.
    assert list(np.asarray(image.get_array())[:, 0]) == [9.0, 2.0, 0.0]
    plt.close("all")


def test_an_all_zero_region_gets_a_unit_range() -> None:
    """The zero branches used to compare bound methods, so they never ran."""
    empty = make_adata(np.zeros((2, 2)), ["chr1", "chr1"], [0, 100])
    plot_region(empty, "chr1")
    track, image = axes_of(plt.gcf())
    assert track.get_ylim() == (0.0, 1.0)
    assert image.get_clim() == (0.0, 1.0)
    plt.close("all")


def test_explicit_limits_win_over_the_data(adata: ad.AnnData) -> None:
    plot_region(adata, "chr1", signal_min=2, signal_max=20, map_min=1, map_max=5)
    track, image = axes_of(plt.gcf())
    assert track.get_ylim() == (2.0, 20.0)
    assert image.get_clim() == (1.0, 5.0)
    plt.close("all")


def test_a_limit_of_zero_is_a_value_not_a_missing_argument(adata: ad.AnnData) -> None:
    """`or` read an explicit 0 as unset and silently used the data bound."""
    plot_region(adata, "chr1", signal_min=0, signal_max=8, map_min=0, map_max=3)
    track, image = axes_of(plt.gcf())
    assert track.get_ylim() == (0.0, 8.0)
    assert image.get_clim() == (0.0, 3.0)
    plt.close("all")


def test_a_dense_matrix_is_accepted_as_well_as_a_sparse_one(adata: ad.AnnData) -> None:
    dense = adata.copy()
    dense.X = np.asarray(adata.X.todense())
    plot_region(dense, "chr1")
    assert axes_of(plt.gcf())[1].get_array().shape == (3, 3)
    plt.close("all")


def test_a_region_with_no_features_is_an_error(adata: ad.AnnData) -> None:
    with pytest.raises(ValueError, match="No features found"):
        plot_region(adata, "chr9")
