"""Tests for ``LSA`` and ``GibbsLDA``, and their use in ``scClusterCells``."""

from __future__ import annotations

from typing import TYPE_CHECKING

import anndata as ad
import numpy as np
import pytest
from _cli_testing import run, tool_path
from scipy import sparse
from scipy.linalg import subspace_angles

from sincei.tools.TopicModels import GENSIM_SLOPE, LSA, GibbsLDA, LSASolver, tfidf

if TYPE_CHECKING:
    from pathlib import Path

N_CELLS = 120
N_REGIONS = 60
N_PC = 5


@pytest.fixture
def counts() -> sparse.csr_matrix:
    """Counts of two groups of cells, each with more counts in half of the regions."""
    rng = np.random.default_rng(0)
    rates = np.full((N_CELLS, N_REGIONS), 0.3)
    rates[: N_CELLS // 2, : N_REGIONS // 2] = 2.0
    rates[N_CELLS // 2 :, N_REGIONS // 2 :] = 2.0
    return sparse.csr_matrix(rng.poisson(rates).astype(np.float32))


def dense_tfidf(matrix: np.ndarray) -> np.ndarray:
    """gensim's lfu TF-IDF, written on the dense matrix."""
    n = matrix.shape[0]
    held = matrix > 0
    holders = held.sum(axis=0)
    region_idf = np.where(holders > 0, np.log2(n / np.maximum(holders, 1)), 0.0)
    weight = np.where(held, 1 + np.log2(np.where(held, matrix, 1)), 0) * region_idf
    pivot = held.sum() / n
    kept = (weight > 0).sum(axis=1, keepdims=True)
    return weight / ((1 - GENSIM_SLOPE) * pivot + GENSIM_SLOPE * kept)


def test_tfidf_is_the_lfu_formula(counts: sparse.csr_matrix) -> None:
    lsa = LSA(N_PC).fit(counts)
    assert lsa.idf_ is not None
    assert lsa.pivot_ is not None

    weighted = tfidf(counts, lsa.idf_, lsa.pivot_).toarray()

    np.testing.assert_allclose(weighted, dense_tfidf(counts.toarray()), rtol=1e-6)


def test_arpack_gives_the_top_singular_vectors(counts: sparse.csr_matrix) -> None:
    lsa = LSA(N_PC, solver="arpack")
    scores = lsa.fit_transform(counts)
    singular = np.linalg.svd(dense_tfidf(counts.toarray()), compute_uv=False)

    assert lsa.solver_ is LSASolver.arpack
    assert scores.shape == (N_CELLS, N_PC)
    assert lsa.components_ is not None
    assert lsa.components_.shape == (N_REGIONS, N_PC)
    np.testing.assert_allclose(lsa.singular_values_, singular[:N_PC], rtol=1e-4)


def test_randomized_and_arpack_find_the_same_components(
    counts: sparse.csr_matrix,
) -> None:
    arpack = LSA(N_PC, solver="arpack").fit(counts)
    randomized = LSA(N_PC, solver="randomized", device="cpu").fit(counts)
    assert arpack.components_ is not None
    assert randomized.components_ is not None

    signal = slice(0, 2)
    angles = subspace_angles(
        arpack.components_[:, signal], randomized.components_[:, signal]
    )
    assert np.degrees(angles).max() < 1.0
    assert randomized.singular_values_ is not None
    assert arpack.singular_values_ is not None
    np.testing.assert_allclose(
        randomized.singular_values_[signal],
        arpack.singular_values_[signal],
        rtol=1e-4,
    )


@pytest.mark.parametrize("solver", ["arpack", "randomized"])
def test_transform_gives_the_scores_of_the_fit(
    counts: sparse.csr_matrix, solver: str
) -> None:
    lsa = LSA(N_PC, solver=solver, device="cpu")
    scores = lsa.fit_transform(counts)

    np.testing.assert_allclose(lsa.transform(counts), scores, rtol=1e-4, atol=1e-5)


def test_the_largest_entry_of_every_component_is_positive(
    counts: sparse.csr_matrix,
) -> None:
    lsa = LSA(N_PC).fit(counts)
    assert lsa.components_ is not None

    largest = np.abs(lsa.components_).argmax(axis=0)
    assert (lsa.components_[largest, np.arange(N_PC)] > 0).all()


def test_lsa_takes_anndata_and_dense_input(counts: sparse.csr_matrix) -> None:
    reference = LSA(N_PC).fit_transform(counts)

    np.testing.assert_allclose(LSA(N_PC).fit_transform(ad.AnnData(counts)), reference)
    np.testing.assert_allclose(LSA(N_PC).fit_transform(counts.toarray()), reference)


def test_lsa_rejects_an_unknown_solver() -> None:
    with pytest.raises(ValueError, match="'arpack'"):
        LSA(N_PC, solver="lanczos")


def test_lsa_rejects_too_many_components(counts: sparse.csr_matrix) -> None:
    with pytest.raises(ValueError, match="n_pc"):
        LSA(N_REGIONS).fit(counts)


def test_lda_gives_topic_proportions(counts: sparse.csr_matrix) -> None:
    lda = GibbsLDA(n_topics=2, iterations=50)
    theta = lda.fit_transform(counts)

    assert theta.shape == (N_CELLS, 2)
    np.testing.assert_allclose(theta.sum(axis=1), 1.0, rtol=1e-5)
    assert lda.topics_ is not None
    assert lda.topics_.shape == (2, N_REGIONS)
    np.testing.assert_allclose(lda.topics_.sum(axis=1), 1.0, rtol=1e-5)
    assert len(lda.log_likelihood_) == 5


def test_lda_separates_the_two_groups(counts: sparse.csr_matrix) -> None:
    theta = GibbsLDA(n_topics=2, iterations=100).fit_transform(counts)

    first, second = theta[: N_CELLS // 2, 0], theta[N_CELLS // 2 :, 0]
    assert abs(first.mean() - second.mean()) > 0.5


def test_lda_is_reproducible_with_a_seed(counts: sparse.csr_matrix) -> None:
    first = GibbsLDA(n_topics=3, iterations=20, random_state=1).fit_transform(counts)
    second = GibbsLDA(n_topics=3, iterations=20, random_state=1).fit_transform(counts)

    np.testing.assert_array_equal(first, second)


def test_a_cell_without_counts_gets_the_prior_mean(counts: sparse.csr_matrix) -> None:
    matrix = counts.tolil()
    matrix[0] = 0
    theta = GibbsLDA(n_topics=4, iterations=10).fit_transform(matrix.tocsr())

    np.testing.assert_allclose(theta[0], 0.25)


def test_lda_rejects_counts_that_are_not_integers(counts: sparse.csr_matrix) -> None:
    with pytest.raises(ValueError, match="integer counts"):
        GibbsLDA(iterations=1).fit(counts * 0.5)

    GibbsLDA(n_topics=2, iterations=1, binarize=True).fit(counts * 0.5)


@pytest.mark.parametrize("method", ["LSA", "LDA"])
def test_scClusterCells_stores_the_reduction_and_the_loadings(
    counts: sparse.csr_matrix, tmp_path: Path, method: str
) -> None:
    source, out = tmp_path / "counts.h5ad", tmp_path / "clustered.h5ad"
    adata = ad.AnnData(counts)
    adata.obs_names = [f"sample_cell{i}" for i in range(N_CELLS)]
    adata.obs["sample"] = "sample"
    adata.obs["barcode"] = [f"cell{i}" for i in range(N_CELLS)]
    adata.var_names = [f"chr1_{i * 1000}_{(i + 1) * 1000}" for i in range(N_REGIONS)]
    adata.var["chrom"] = "chr1"
    adata.var["name"] = adata.var_names
    adata.var["start"] = np.arange(N_REGIONS) * 1000
    adata.var["end"] = adata.var["start"] + 1000
    adata.write_h5ad(source)

    proc = run(
        tool_path("scClusterCells"),
        ["-i", str(source), "-o", str(out), "-m", method, "-n", "3", "-nk", "10"]
        + (["--nIterations", "20"] if method == "LDA" else []),
    )
    assert proc.returncode == 0, proc.stderr

    result = ad.read_h5ad(out)
    assert result.obsm[f"X_{method}"].shape == (N_CELLS, 3)
    assert result.varm[f"X_{method}"].shape == (N_REGIONS, 3)
    assert "leiden" in result.obs
