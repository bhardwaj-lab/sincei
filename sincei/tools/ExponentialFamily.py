from __future__ import annotations

import copy
from enum import Enum
from typing import Any

import numpy as np
import torch
import torch.distributions as td

EULER_GAMMA = 0.5772156649015329


class GLMFamily(str, Enum):
    gaussian = "gaussian"
    poisson = "poisson"
    negative_binomial = "negative_binomial"
    bernoulli = "bernoulli"
    beta = "beta"
    gamma = "gamma"
    lognormal = "lognormal"
    sigmoid_beta = "sigmoid_beta"

    def distribution(self) -> type[ExponentialFamily]:
        return {
            GLMFamily.gaussian: Gaussian,
            GLMFamily.poisson: Poisson,
            GLMFamily.negative_binomial: NegativeBinomial,
            GLMFamily.bernoulli: Bernoulli,
            GLMFamily.beta: Beta,
            GLMFamily.gamma: Gamma,
            GLMFamily.lognormal: LogNormal,
            GLMFamily.sigmoid_beta: SigmoidBeta,
        }[self]


def inverse_digamma(y: torch.Tensor, n_iter: int = 8) -> torch.Tensor:
    """Vectorised inverse digamma (Minka, 2000): asymptotic init + Newton steps."""
    x = torch.where(y >= -2.22, torch.exp(y) + 0.5, -1.0 / (y + EULER_GAMMA))
    for _ in range(n_iter):
        x = (x - (torch.digamma(x) - y) / torch.polygamma(1, x)).clamp_min(1e-12)
    return x


class ExponentialFamily:
    r"""Exponential family distribution for GLM-PCA.

    A subclass needs two things:

    1. The likelihood, given in one of two ways:
       - ``torch_distribution(theta)``: return a ``torch.distributions`` object
         parametrised by theta (preferred: exact log-likelihood incl. base measure,
         nothing to derive), or
       - the paper's notation (Mourragui et al., 2023): ``sufficient_statistics``
         (T), ``natural_parametrization`` (eta) and ``log_partition`` (A).
    2. ``invert_g``: the element-wise saturated parameter g^{-1}(x).

    Optional: per-feature nuisance parameters (e.g. dispersion) listed in
    ``nuisance_params`` and estimated from streamed column sums
    (``family_statistics`` / ``initialize_family_parameters_from_statistics``).

    Convention: any 1-D tensor in ``family_params`` holds one value per feature,
    which lets the driver subset features generically (``subset_features``).
    """

    family_name: str = "base"
    default_params: dict[str, Any] = {}
    nuisance_params: tuple[str, ...] = ()

    def __init__(
        self, family_params: dict[str, Any] | None = None, **kwargs: object
    ) -> None:
        self.family_params = {**self.default_params, **(family_params or {}), **kwargs}

    # ---- likelihood -------------------------------------------------------------
    def torch_distribution(self, theta: torch.Tensor) -> td.Distribution:
        raise NotImplementedError

    def sufficient_statistics(self, X: torch.Tensor) -> torch.Tensor:
        return X

    def natural_parametrization(self, theta: torch.Tensor) -> torch.Tensor:
        return theta

    def log_partition(self, theta: torch.Tensor) -> torch.Tensor:
        raise NotImplementedError

    def clip_data(self, X: torch.Tensor) -> torch.Tensor:
        """Maps data into the support of the distribution (identity by default)."""
        return X

    def invert_g(self, X: torch.Tensor) -> torch.Tensor:
        return X

    def log_likelihood(self, X: torch.Tensor, theta: torch.Tensor) -> torch.Tensor:
        """Negative log-likelihood summed over all entries."""
        X = self.clip_data(X)
        if type(self).torch_distribution is not ExponentialFamily.torch_distribution:
            return -self.torch_distribution(theta).log_prob(X).sum()
        T, eta = self.sufficient_statistics(X), self.natural_parametrization(theta)
        expt = (T * eta).sum(0) if T.dim() > X.dim() else T * eta
        return -(expt - self.log_partition(theta)).sum()

    # ---- per-feature parameters -------------------------------------------------
    def _param(self, name: str, like: torch.Tensor) -> torch.Tensor:
        return torch.as_tensor(self.family_params[name]).to(like)

    def to(self, device: torch.device | str, dtype: torch.dtype | None = None):
        self.family_params = {
            k: v.to(device, dtype if v.is_floating_point() else None)
            if isinstance(v, torch.Tensor) else v
            for k, v in self.family_params.items()
        }
        return self

    def load_family_params_to_gpu(self, device: torch.device) -> None:
        self.to(device)

    def subset_features(self, idx: torch.Tensor) -> ExponentialFamily:
        """Shallow copy with every 1-D (per-feature) parameter indexed by ``idx``
        (a feature subset, or one feature index per non-zero entry)."""
        new = copy.copy(self)
        new.family_params = {
            k: v[idx] if isinstance(v, torch.Tensor) and v.dim() == 1 else v
            for k, v in self.family_params.items()
        }
        return new

    # ---- streaming nuisance-parameter estimation --------------------------------
    def needs_family_statistics(self) -> bool:
        return any(k not in self.family_params for k in self.nuisance_params)

    def family_statistics(self, X: torch.Tensor) -> dict[str, torch.Tensor]:
        """Per-feature column sums of a chunk (summed over chunks by the caller)."""
        return {}

    def initialize_family_parameters_from_statistics(
        self, stats: dict[str, torch.Tensor], n: int
    ) -> None:
        """Sets nuisance parameters from accumulated column sums over ``n`` rows."""

    def initialize_family_parameters(self, X: torch.Tensor) -> None:
        if self.needs_family_statistics():
            self.initialize_family_parameters_from_statistics(
                self.family_statistics(X), X.shape[0]
            )


# ======================================================================================
class Gaussian(ExponentialFamily):
    """Unit-variance Gaussian: GLM-PCA reduces to PCA."""

    family_name = "gaussian"

    def torch_distribution(self, theta):
        return td.Normal(theta, 1.0, validate_args=False)


class Bernoulli(ExponentialFamily):
    """Bernoulli. ``max_val``: truncation m of the infinite saturated parameters
    (Landgraf & Lee, 2020); acts as a regulariser, tune by held-out deviance."""

    family_name = "bernoulli"
    default_params = {"max_val": 30}

    def torch_distribution(self, theta):
        return td.Bernoulli(logits=theta, validate_args=False)

    def invert_g(self, X):
        m = self.family_params["max_val"]
        eps = 1.0 / (1.0 + np.exp(m))
        return torch.logit(X.clamp(eps, 1 - eps)).clamp(-m, m)


class Poisson(ExponentialFamily):
    """Poisson, theta = log(rate). ``min_val`` replaces 0 before the log."""

    family_name = "poisson"
    default_params = {"min_val": 1e-20}

    def torch_distribution(self, theta):
        return td.Poisson(rate=theta.exp(), validate_args=False)

    def invert_g(self, X):
        return torch.log(X.clamp_min(self.family_params["min_val"]))


class NegativeBinomial(Poisson):
    """Negative binomial, theta = log(mean), fixed per-feature inverse dispersion r
    (Risso et al.). r is estimated by streaming method of moments if not given."""

    family_name = "negative_binomial"
    default_params = {"min_val": 1e-20, "r_min": 1e-2, "r_max": 1e4}
    nuisance_params = ("r",)

    def torch_distribution(self, theta):
        r = self._param("r", theta)
        # torch's NB: mean = r * exp(logits)  ->  logits = theta - log r
        return td.NegativeBinomial(r, logits=theta - r.log(), validate_args=False)

    def family_statistics(self, X):
        X = X.double()
        return {"s1": X.sum(0), "s2": (X * X).sum(0)}

    def initialize_family_parameters_from_statistics(self, stats, n):
        m = stats["s1"] / n
        v = stats["s2"] / n - m**2
        r_max = self.family_params["r_max"]
        r = torch.where(v > m, m**2 / (v - m).clamp_min(1e-12), torch.full_like(m, r_max))
        self.family_params["r"] = r.clamp(self.family_params["r_min"], r_max).float()


class Beta(ExponentialFamily):
    """Beta with mean parametrisation theta in (0, 1) and fixed per-feature
    concentration nu (estimated by streaming MLE if not given)."""

    family_name = "beta"
    default_params = {"min_val": 1e-5, "eps": 1e-8, "maxiter": 100, "method": "MLE"}
    nuisance_params = ("nu",)

    def clip_data(self, X):
        mv = self.family_params["min_val"]
        return X.clamp(mv, 1 - mv)

    def _mean(self, theta):
        # the linear projection can leave (0, 1); prefer SigmoidBeta
        return self.clip_data(theta)

    def torch_distribution(self, theta):
        m, nu = self._mean(theta), self._param("nu", theta)
        return td.Beta(m * nu, (1 - m) * nu, validate_args=False)

    def _saturated_mean(self, X):
        """Solves logit(x) + psi((1-m) nu) - psi(m nu) = 0 (decreasing in m)."""
        from scipy.optimize.elementwise import find_root
        from scipy.special import digamma, logit

        x = self.clip_data(X).double().cpu().numpy()
        nu = self._param("nu", X).double().cpu().numpy()
        res = find_root(
            lambda m, lx, nu: lx + digamma((1 - m) * nu) - digamma(m * nu),
            (np.full_like(x, 1e-12), np.full_like(x, 1 - 1e-12)),
            args=(logit(x), nu),
            tolerances={"xatol": self.family_params["eps"]},
            maxiter=self.family_params["maxiter"],
        )
        return torch.as_tensor(res.x).to(X)

    def invert_g(self, X):
        return self._saturated_mean(X)

    def family_statistics(self, X):
        X = self.clip_data(X).double()
        return {"s1": X.sum(0), "s2": (X * X).sum(0),
                "slx": X.log().sum(0), "sl1x": torch.log1p(-X).sum(0)}

    def initialize_family_parameters_from_statistics(self, stats, n):
        m = stats["s1"] / n
        v = (stats["s2"] / n - m**2).clamp_min(1e-12)
        c = (m * (1 - m) / v - 1).clamp_min(1e-2)  # method of moments
        a, b = (m * c).clamp_min(1e-3), ((1 - m) * c).clamp_min(1e-3)
        if self.family_params["method"] == "MLE":  # Minka fixed point, all features
            mlx, ml1x = stats["slx"] / n, stats["sl1x"] / n
            for _ in range(200):
                s = torch.digamma(a + b)
                a_new, b_new = inverse_digamma(s + mlx), inverse_digamma(s + ml1x)
                done = max((a_new - a).abs().max(), (b_new - b).abs().max()) < 1e-8
                a, b = a_new, b_new
                if done:
                    break
        self.family_params["nu"] = (a + b).float()


class SigmoidBeta(Beta):
    """Beta with mean = sigmoid(theta): unconstrained theta. Note this changes the
    model (linear structure is sought in logit(mean))."""

    family_name = "sigmoid_beta"

    def _mean(self, theta):
        return torch.sigmoid(theta)

    def invert_g(self, X):
        return torch.logit(self._saturated_mean(X))


class Gamma(ExponentialFamily):
    """Gamma with shape theta + 1 and fixed per-feature rate nu (as in the paper).
    Saturated theta = psi^{-1}(log(nu x)) - 1."""

    family_name = "gamma"
    default_params = {"min_val": 1e-5, "max_val": 1e7}
    nuisance_params = ("nu",)

    def clip_data(self, X):
        return X.clamp(self.family_params["min_val"], self.family_params["max_val"])

    def torch_distribution(self, theta):
        return td.Gamma(theta + 1, self._param("nu", theta), validate_args=False)

    def invert_g(self, X):
        return inverse_digamma(torch.log(self.clip_data(X) * self._param("nu", X))) - 1.0

    def family_statistics(self, X):
        X = self.clip_data(X).double()
        return {"s1": X.sum(0), "slx": X.log().sum(0)}

    def initialize_family_parameters_from_statistics(self, stats, n):
        mean = stats["s1"] / n
        s = (torch.log(mean) - stats["slx"] / n).clamp_min(1e-10)
        a = (3 - s + torch.sqrt((s - 3) ** 2 + 24 * s)) / (12 * s)  # Minka init
        for _ in range(20):  # Newton on log a - psi(a) = s
            a = (a - (a.log() - torch.digamma(a) - s)
                 / (1 / a - torch.polygamma(1, a))).clamp_min(1e-6)
        self.family_params["nu"] = (a / mean).float()  # rate


class LogNormal(ExponentialFamily):
    """Log-normal with fixed per-feature sigma nu. Saturated theta = log x."""

    family_name = "lognormal"
    default_params = {"min_val": 1e-5}
    nuisance_params = ("nu",)

    def clip_data(self, X):
        return X.clamp_min(self.family_params["min_val"])

    def torch_distribution(self, theta):
        return td.LogNormal(theta, self._param("nu", theta), validate_args=False)

    def invert_g(self, X):
        return torch.log(self.clip_data(X))

    def family_statistics(self, X):
        lx = torch.log(self.clip_data(X)).double()
        return {"s1": lx.sum(0), "s2": (lx * lx).sum(0)}

    def initialize_family_parameters_from_statistics(self, stats, n):
        var = (stats["s2"] - stats["s1"] ** 2 / n) / max(n - 1, 1)
        self.family_params["nu"] = var.clamp_min(1e-8).sqrt().float()
