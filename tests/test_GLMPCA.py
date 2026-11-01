"""Tests for ``GLMPCA``: family selection, the fit/transform contract, AnnData input."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import anndata as ad
import numpy as np
import pytest
import torch
from scipy import sparse

from sincei.tools.ExponentialFamily import GLMFamily, Poisson
from sincei.tools.GLMPCA import GLMPCA

if TYPE_CHECKING:
    from collections.abc import Callable

N_CELLS = 40
N_FEATURES = 12
N_PC = 2


@pytest.fixture(autouse=True)
def seed() -> None:
    torch.manual_seed(0)
    np.random.seed(0)


def sample(family: GLMFamily) -> torch.Tensor:
    """Data inside the support of ``family``."""
    shape = (N_CELLS, N_FEATURES)
    if family is GLMFamily.gaussian:
        return torch.randn(shape)
    if family is GLMFamily.poisson:
        return torch.poisson(torch.full(shape, 3.0))
    if family is GLMFamily.bernoulli:
        return torch.bernoulli(torch.full(shape, 0.4))
    if family in {GLMFamily.beta, GLMFamily.sigmoid_beta}:
        return torch.rand(shape) * 0.8 + 0.1
    return torch.rand(shape) + 0.5


@pytest.mark.parametrize("family", list(GLMFamily))
def test_fit_and_transform_project_onto_n_pc_components(family: GLMFamily) -> None:
    X = sample(family)
    model = GLMPCA(N_PC, family=family, max_iter=2, batch_size=16)

    assert model.fit(X)
    assert model.saturated_loadings_ is not None
    assert model.saturated_loadings_.shape == (N_FEATURES, N_PC)
    assert model.transform(X).shape == (N_CELLS, N_PC)


@pytest.mark.parametrize("family", list(GLMFamily))
def test_a_family_name_and_its_member_select_the_same_distribution(
    family: GLMFamily,
) -> None:
    assert type(GLMPCA(N_PC, family=family.value).exponential_family) is (
        family.distribution()
    )
    assert type(GLMPCA(N_PC, family=family).exponential_family) is (
        family.distribution()
    )


def test_an_exponential_family_instance_is_used_as_given() -> None:
    family = Poisson({"min_val": 1e-10})
    assert GLMPCA(N_PC, family=family).exponential_family is family


def test_an_unknown_family_name_is_rejected() -> None:
    with pytest.raises(ValueError, match="nope"):
        GLMPCA(N_PC, family="nope")


def test_transform_before_fit_fails_with_advice() -> None:
    with pytest.raises(RuntimeError, match=r"Call fit\(\)"):
        GLMPCA(N_PC, family="poisson").transform(sample(GLMFamily.poisson))


def fitted_loadings(X: torch.Tensor | ad.AnnData) -> torch.Tensor:
    model = GLMPCA(N_PC, family="poisson", max_iter=2, batch_size=8, seed=0)
    model.fit(X)
    assert model.saturated_loadings_ is not None
    return model.saturated_loadings_.detach()


@pytest.mark.parametrize(
    "matrix_type",
    [
        np.asarray,
        sparse.csr_matrix,
        sparse.csc_matrix,
        sparse.csr_array,
        sparse.csc_array,
    ],
    ids=lambda matrix_type: matrix_type.__name__,
)
def test_anndata_input_is_fitted_with_cells_in_rows(
    matrix_type: Callable[[np.ndarray], Any],
) -> None:
    counts = (
        np.random
        .default_rng(0)
        .poisson(3.0, size=(N_CELLS, N_FEATURES))
        .astype(np.float32)
    )
    adata = ad.AnnData(matrix_type(counts))

    # sparse input takes the theta0 + Delta path, dense the direct one: same model
    torch.testing.assert_close(
        fitted_loadings(adata), fitted_loadings(torch.Tensor(counts)),
        rtol=1e-4, atol=1e-4,
    )


@pytest.mark.parametrize("family", list(GLMFamily))
def test_sparse_and_dense_input_give_the_same_scores(family: GLMFamily) -> None:
    X = sample(family) * (torch.rand(N_CELLS, N_FEATURES) > 0.5)
    model = GLMPCA(N_PC, family=family, max_iter=2, batch_size=16, seed=0)
    model.fit(X)

    torch.testing.assert_close(
        model.transform(sparse.csr_matrix(X.numpy())), model.transform(X),
        rtol=1e-4, atol=1e-4,
    )


def test_refit_re_estimates_per_feature_parameters() -> None:
    model = GLMPCA(N_PC, family="negative_binomial", max_iter=2, seed=0)
    model.fit(sample(GLMFamily.poisson))
    model.fit(sample(GLMFamily.poisson)[:, :-2])  # fewer features

    assert model.exponential_family.family_params["r"].shape == (N_FEATURES - 2,)


def _grads(model: GLMPCA, X: Any, n_jobs: int = 1) -> tuple[float, torch.Tensor, torch.Tensor]:
    """Cost and gradients of one full pass, at the fitted point."""
    model.n_jobs = n_jobs
    data = model._data(X)
    V = model.saturated_loadings_.clone().requires_grad_(True)
    mu = model.saturated_intercept_.clone().requires_grad_(True)
    cost = model._accumulate(V, mu, data, np.arange(data.n), backward=True, scale=1 / data.n)
    return cost, V.grad, mu.grad


def test_sparse_path_gives_the_dense_gradient_with_an_offset() -> None:
    X = sample(GLMFamily.poisson) * (torch.rand(N_CELLS, N_FEATURES) > 0.5)
    model = GLMPCA(N_PC, family="poisson", max_iter=2, seed=0, depth_correction="offset")
    model.fit(X)
    Xs = sparse.csr_matrix(X.numpy())
    assert model._data(Xs).batch(np.arange(4)).delta.is_coalesced()

    for got, want in zip(_grads(model, Xs), _grads(model, X)):
        torch.testing.assert_close(got, want, rtol=1e-4, atol=1e-4)


def test_threaded_passes_match_the_serial_pass() -> None:
    X = sparse.csr_matrix(sample(GLMFamily.poisson).numpy())
    model = GLMPCA(N_PC, family="poisson", max_iter=2, seed=0, chunk_size=8)  # 4 pieces of 2
    model.fit(X)

    for got, want in zip(_grads(model, X, n_jobs=4), _grads(model, X, n_jobs=1)):
        torch.testing.assert_close(got, want, rtol=1e-5, atol=1e-5)


@pytest.mark.parametrize("optimizer", ["geoopt", "lbfgs"])
def test_n_jobs_fits_and_restores_torch_threads(optimizer: str) -> None:
    threads = torch.get_num_threads()
    X = sample(GLMFamily.poisson)
    fits = [GLMPCA(N_PC, family="poisson", optimizer=optimizer, fit_mode="exact",
                   max_iter=3, chunk_size=8, n_jobs=n_jobs, seed=0).fit_transform(X)
            for n_jobs in (1, 3)]

    torch.testing.assert_close(fits[1], fits[0], rtol=1e-3, atol=1e-3)
    assert torch.get_num_threads() == threads


def test_prefetching_minibatches_does_not_change_the_fit() -> None:
    X = sample(GLMFamily.poisson)
    fits = []
    for prefetch in (True, False):
        model = GLMPCA(N_PC, family="poisson", max_iter=3, batch_size=8,
                       feature_subsample=0.5, seed=0)
        model._prefetch = prefetch
        model.fit(X)
        fits.append(model.saturated_loadings_)

    assert torch.equal(fits[0], fits[1])


@pytest.mark.parametrize("optimizer", ["geoopt", "lbfgs"])
def test_threaded_fits_with_the_same_seed_are_identical(optimizer: str) -> None:
    X = sample(GLMFamily.poisson)
    fits = [GLMPCA(N_PC, family="poisson", optimizer=optimizer, fit_mode="exact",
                   max_iter=3, chunk_size=8, n_jobs=4, seed=0).fit_transform(X)
            for _ in range(2)]

    assert torch.equal(fits[0], fits[1])
