"""Check each exponential family against the scipy distribution it encodes.

Each test sweeps a grid of parameters and compares the density, or its log, over
a grid of observations.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
import scipy
import torch

from sincei.tools.ExponentialFamily import (
    Bernoulli,
    Beta,
    Gamma,
    Gaussian,
    NegativeBinomial,
    Poisson,
    SigmoidBeta,
)

if TYPE_CHECKING:
    from collections.abc import Callable


@pytest.fixture(autouse=True)
def seed() -> None:
    torch.manual_seed(0)


def test_gaussian_pdf_matches_scipy() -> None:
    X = torch.linspace(-5, 5, 500)
    for theta in torch.normal(0.0, 1.0, (25,)):
        pdf = Gaussian().torch_distribution(theta.expand_as(X)).log_prob(X).exp()
        np.testing.assert_array_almost_equal(
            pdf.numpy(), scipy.stats.norm.pdf(X.numpy(), loc=float(theta)), decimal=3
        )


def test_bernoulli_success_probability_is_the_sigmoid_of_theta() -> None:
    theta = torch.linspace(-30, 30, 1000)
    success = Bernoulli().torch_distribution(theta).log_prob(torch.ones_like(theta)).exp()
    torch.testing.assert_close(success, torch.sigmoid(theta), rtol=0, atol=1e-6)


def test_poisson_log_pmf_matches_scipy() -> None:
    X = torch.linspace(0, 100, 101)
    for theta in torch.linspace(-50, 5, 10):
        log_pmf = Poisson().torch_distribution(theta.expand_as(X)).log_prob(X)
        np.testing.assert_array_almost_equal(
            log_pmf.numpy(),
            scipy.stats.poisson.logpmf(X.numpy(), np.exp(float(theta))),
            decimal=3,
        )


@pytest.mark.parametrize(
    ("family", "thetas", "mean"),
    [
        (Beta, torch.linspace(0.01, 0.99, 10), lambda theta: theta),
        (SigmoidBeta, torch.logit(torch.linspace(0.01, 0.99, 10)), torch.sigmoid),
    ],
    ids=["beta", "sigmoid_beta"],
)
def test_beta_log_pdf_matches_scipy(
    family: type[Beta],
    thetas: torch.Tensor,
    mean: Callable[[torch.Tensor], torch.Tensor],
) -> None:
    X = torch.linspace(1e-4, 1 - 1e-4, 100)
    for nu in torch.rand(20) * 10:
        for theta in thetas:
            distribution = family()
            distribution.family_params["nu"] = nu.expand_as(X)
            log_pdf = distribution.torch_distribution(theta.expand_as(X)).log_prob(X)
            p = float(mean(theta))
            np.testing.assert_array_almost_equal(
                log_pdf.numpy(),
                scipy.stats.beta.logpdf(X.numpy(), p * float(nu), (1 - p) * float(nu)),
                decimal=2,
            )


def test_gamma_pdf_matches_scipy() -> None:
    X = torch.logspace(-8, 5, 1000)
    for nu in torch.rand(20) * 10:
        for theta in torch.logspace(-5, 3, 50):
            distribution = Gamma()
            distribution.family_params["nu"] = nu.expand_as(X)
            pdf = distribution.torch_distribution(theta.expand_as(X)).log_prob(X).exp()
            np.testing.assert_array_almost_equal(
                pdf.numpy(),
                scipy.stats.gamma.pdf(
                    X.numpy(), a=float(theta) + 1, loc=0, scale=1 / float(nu)
                ),
                decimal=2,
            )


@pytest.mark.parametrize(
    "family", [Poisson(), NegativeBinomial({"r": torch.tensor([0.5, 2.0, 50.0])})],
    ids=lambda family: family.family_name,
)
def test_fit_nll_differs_from_the_full_nll_by_a_constant(family) -> None:
    X = torch.poisson(torch.full((30, 3), 4.0))
    gaps = []
    for _ in range(2):
        theta = torch.randn(30, 3, requires_grad=True)
        full, fit = family.log_likelihood(X, theta), family.fit_nll(X, theta)
        g_full, = torch.autograd.grad(full, theta)
        g_fit, = torch.autograd.grad(fit, theta)
        torch.testing.assert_close(g_fit, g_full)
        gaps.append(float(full - fit))
    assert gaps[0] == pytest.approx(gaps[1], rel=1e-5)
