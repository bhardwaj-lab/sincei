import torch
from copy import deepcopy
import numpy as np
import scipy
from tqdm import tqdm
from joblib import Parallel, delayed

saturation_eps = 10**-10


class ExponentialFamily:
    r"""Encodes an exponential family distribution using PyTorch autodiff structures.

    ExponentialFamily corresponds to the superclass providing a backbone for
    all exponential family.
    Each subclass should contain the following methods, defined based on the
    distribution of choice (same notation as [Mourragui et al, 2023]):
        - sufficient_statistics ($T$),
        - natural_parametrization ($\eta$),
        - log_partition ($A$).
        - invert_g ($g^{-1}$).
        - initialize_family_parameters: computes parameters used in other
        methods, e.g., gene-level dispersion for Negative Binomial.
    We added a "base_measure" for sake of completeness, but this method is not
    necessary for running GLM-PCA.
    The log-likelihood and exponential term are defined directly from the
    aforementionned methods.

    Parameters
    ----------
    family_name : int
        Name of the family.
    """

    def __init__(self, family_params=None, **kwargs):
        self.family_name = "base"
        self.family_params = family_params if family_params else {}

    def sufficient_statistics(self, X: torch.Tensor):
        return X

    def natural_parametrization(self, theta: torch.Tensor):
        return theta

    def log_partition(self, theta: torch.Tensor):
        return 0.0

    def base_measure(self, X: torch.Tensor):
        return torch.ones(size=X.shape)

    def invert_g(self, X: torch.Tensor):
        return X

    def exponential_term(self, X: torch.Tensor, theta: torch.Tensor):
        return torch.multiply(self.sufficient_statistics(X), self.natural_parametrization(theta))

    def distribution(self, X: torch.Tensor, theta: torch.Tensor):
        f = self.base_measure(X)
        expt = self.exponential_term(X, theta) - self.log_partition(theta)

        return torch.multiply(f, torch.exp(expt))

    def log_distribution(self, X: torch.Tensor, theta: torch.Tensor):
        f = self.base_measure(X)
        expt = self.exponential_term(X, theta) - self.log_partition(theta)

        return expt - torch.log(f)

    def log_likelihood(self, X: torch.Tensor, theta: torch.Tensor):
        """Computes negative log-likelihood between dataset X and parameters theta"""
        expt = self.exponential_term(X, theta) - self.log_partition(theta)
        return -torch.sum(expt)

    def load_family_params_to_gpu(self, device):
        self.family_params = {
            k: (
                self.family_params[k].to(device)
                if type(self.family_params[k]) is torch.Tensor
                else self.family_params[k]
            )
            for k in self.family_params
        }

    def initialize_family_parameters(self, X: torch.Tensor = None):
        """General method to initialize certain parameters (e.g. for Beta or Negative Binomial"""
        pass


class Gaussian(ExponentialFamily):
    r"""Gaussian with standard deviation one.

    GLMPCA with Gaussian as family is equivalent to the standard PCA.
    """

    def __init__(self, family_params=None, **kwargs):
        self.family_name = "gaussian"
        self.family_params = family_params if family_params else {}

    def sufficient_statistics(self, X: torch.Tensor):
        return X

    def natural_parametrization(self, theta: torch.Tensor):
        return theta

    def log_partition(self, theta: torch.Tensor):
        return torch.square(theta) / 2.0

    def base_measure(self, X: torch.Tensor):
        return torch.exp(-torch.square(X) / 2.0) / np.sqrt(2.0 * torch.pi)

    def invert_g(self, X: torch.Tensor):
        return X


class Bernoulli(ExponentialFamily):
    r"""Bernoulli distribution

    family_params of interest:
        - "max_val" (int) corresponding to the max value (replaces infinity).
        Empirically, values above 10 yield similar results.
    """

    def __init__(self, family_params=None, **kwargs):
        self.family_name = "bernoulli"
        self.family_params = family_params if family_params else {"max_val": 30}

    def sufficient_statistics(self, X: torch.Tensor):
        return X

    def natural_parametrization(self, theta: torch.Tensor):
        return theta

    def log_partition(self, theta: torch.Tensor):
        return torch.log(1.0 + torch.exp(theta))

    def base_measure(self, X: torch.Tensor):
        return 1.0

    def invert_g(self, X: torch.Tensor):
        return torch.log(X / (X - 1)).clip(-self.family_params["max_val"], self.family_params["max_val"])

    def log_likelihood(self, X: torch.Tensor, theta: torch.Tensor):
        """Computes negative log-likelihood between dataset X and parameters theta"""
        expt = self.exponential_term(X, theta) - self.log_partition(theta)
        return -torch.sum(expt)


class Poisson(ExponentialFamily):
    r"""Poisson distribution

    family_params of interest:
        - "min_val" (int) corresponding to the min value (replaces 0).
    """

    def __init__(self, family_params=None, **kwargs):
        self.family_name = "poisson"
        self.family_params = family_params if family_params else {"min_val": 1e-20}

    def sufficient_statistics(self, X: torch.Tensor):
        return X

    def natural_parametrization(self, theta: torch.Tensor):
        return theta

    def log_partition(self, theta: torch.Tensor):
        return torch.exp(theta)

    def base_measure(self, X: torch.Tensor):
        return 1 / scipy.special.gamma(X + 1)

    def invert_g(self, X: torch.Tensor):
        return torch.log(X).clip(self.family_params["min_val"])

    def log_distribution(self, X: torch.Tensor, theta: torch.Tensor):
        """The computation of gamma function for the base measure (h) would lead to inf, hence a re-design
        of the method."""
        log_f = torch.lgamma(X + 1)
        expt = self.exponential_term(X, theta) - self.log_partition(theta)

        return expt - log_f


class Beta(ExponentialFamily):
    r"""Beta distribution, using a standard formulation.

    Original formulation presented in [Mourragui et al, 2023].

    family_params of interest:
        - "min_val" (int): min data value (replaces 0 and 1).
        - "n_jobs" (int): number of jobs, specifically
        for computing the "nu" parameter.
        - "method" (str): method use to compute the "nu" parameter per
        feature. Two possibles: "MLE" and "MM". Default to "MLE".
        - "eps" (float): minimum difference used for inverting the g
        function. Default to 1e-4
        - "maxiter" (int): maximum number of iterations for the inversion
        of the g function. Default to 100.
    """

    def __init__(self, family_params=None, **kwargs):
        self.family_name = "beta"
        if family_params is None or "nu" not in family_params:
            print("Beta distribution not initialized yet")
        default_family_params = {"min_val": 1e-5, "n_jobs": 1, "eps": 1e-4, "maxiter": 100, "method": "MLE"}
        self.family_params = family_params if family_params else default_family_params
        for k in kwargs:
            self.family_params[k] = kwargs[k]
        for k in default_family_params:
            if k in self.family_params:
                continue
            self.family_params[k] = default_family_params[k]

    def sufficient_statistics(self, X: torch.Tensor):
        X = X.clip(self.family_params["min_val"], 1 - self.family_params["min_val"])
        return torch.stack([torch.log(X), torch.log(1 - X)])

    def natural_parametrization(self, theta: torch.Tensor):
        nat_params = torch.stack([theta * self.family_params["nu"], (1 - theta) * self.family_params["nu"]])
        if nat_params.shape[1] == 1:
            nat_params = nat_params.flatten()
        return nat_params

    def base_measure(self, X: torch.Tensor):
        return torch.mul(X, 1 - X)

    def log_partition(self, theta: torch.Tensor):
        numerator = torch.sum(torch.lgamma(self.natural_parametrization(theta)), axis=0)
        # numerator = torch.sum(numerator, axis=0)
        denominator = torch.lgamma(self.family_params["nu"])
        return numerator - denominator

    def exponential_term(self, X: torch.Tensor, theta: torch.Tensor):
        return torch.sum(torch.multiply(self.sufficient_statistics(X), self.natural_parametrization(theta)), axis=0)

    def _derivative_log_likelihood(self, X: torch.Tensor, theta: torch.Tensor):
        X = X.clip(self.family_params["eps"], 1 - self.family_params["eps"])
        return (
            torch.log(X / (1 - X))
            + torch.digamma((1 - theta) * self.family_params["nu"])
            - torch.digamma(theta * self.family_params["nu"])
        )

    def initialize_family_parameters(self, X: torch.Tensor = None):
        p = X.shape[1]

        def compute_beta_param(x):
            y = x[x > self.family_params["eps"]]
            y = y[y < 1 - self.family_params["eps"]]
            return scipy.stats.beta.fit(y, floc=0, fscale=1, method=self.family_params["method"])

        self.family_params["nu"] = torch.Tensor(
            Parallel(n_jobs=self.family_params["n_jobs"], batch_size=100, backend="threading")(
                delayed(compute_beta_param)(X[:, idx]) for idx in tqdm(range(p))
            )
        )
        self.family_params["nu"] = torch.sum(self.family_params["nu"][:, :2], axis=1)
        assert self.family_params["nu"].shape[0] == p

        return True

    def invert_g(self, X: torch.Tensor = None):
        """Dichotomy to find where derivative maxes out"""

        X = X.clip(self.family_params["eps"], 1 - self.family_params["eps"])

        # Initialize dichotomy parameters.
        min_val = torch.zeros(X.shape)
        max_val = torch.ones(X.shape)
        theta = (min_val + max_val) / 2

        llik = self._derivative_log_likelihood(X, theta)
        for idx in tqdm(range(self.family_params["maxiter"])):
            min_val[llik > 0] = theta[llik > 0]
            max_val[llik < 0] = theta[llik < 0]
            theta = (min_val + max_val) / 2
            llik = self._derivative_log_likelihood(X, theta)

            if torch.max(torch.abs(llik)) < self.family_params["eps"]:
                print("CONVERGENCE AFTER %s ITERATIONS" % (idx))
                break

        if idx <= self.family_params["maxiter"]:
            print("CONVERGENCE NOT REACHED")

        return theta


class SigmoidBeta(Beta):
    r"""Beta distribution re-parametrized using a Sigmoid.

    This distribution is similar to the previous Beta (which it
    inherits from) but the natural parameter is re-parametrized using
    a Sigmoid. This is shown expeerimentally to stabilize the
    optimisation by removing the ]0,1[ constraint.

    family_params of interest:
        - "min_val" (int): min data value (replaces 0 and 1).
        - "n_jobs" (int): number of jobs, specifically
        for computing the "nu" parameter.
        - "method" (str): method use to compute the "nu" parameter per
        feature. Two possibles: "MLE" and "MM". Default to "MLE".
        - "eps" (float): minimum difference used for inverting the g
        function. Default to 1e-4
        - "maxiter" (int): maximum number of iterations for the inversion
        of the g function. Default to 100.
    """

    def natural_parametrization(self, theta: torch.Tensor):
        nat_params = torch.stack(
            [torch.sigmoid(theta) * self.family_params["nu"], (1 - torch.sigmoid(theta)) * self.family_params["nu"]]
        )
        if nat_params.shape[1] == 1:
            nat_params = nat_params.flatten()
        return nat_params

    def _derivative_log_likelihood(self, X: torch.Tensor, logit_theta: torch.Tensor):
        X = X.clip(self.family_params["eps"], 1 - self.family_params["eps"])
        return (
            torch.log(X / (1 - X))
            + torch.digamma((1 - logit_theta) * self.family_params["nu"])
            - torch.digamma(logit_theta * self.family_params["nu"])
        )

    def invert_g(self, X: torch.Tensor = None):
        """Dichotomy to find where derivative maxes out"""

        X = X.clip(self.family_params["eps"], 1 - self.family_params["eps"])

        # Initialize dichotomy parameters.
        min_val = torch.zeros(X.shape)
        max_val = torch.ones(X.shape)
        logit_theta = (min_val + max_val) / 2

        llik = self._derivative_log_likelihood(X, logit_theta)
        for idx in tqdm(range(self.family_params["maxiter"])):
            min_val[llik > 0] = logit_theta[llik > 0]
            max_val[llik < 0] = logit_theta[llik < 0]
            logit_theta = (min_val + max_val) / 2
            llik = self._derivative_log_likelihood(X, logit_theta)

            if torch.max(torch.abs(llik)) < self.family_params["eps"]:
                print("CONVERGENCE AFTER %s ITERATIONS" % (idx))
                break

        if idx <= self.family_params["maxiter"]:
            print("CONVERGENCE NOT REACHED")

        return torch.logit(logit_theta)


class Gamma(ExponentialFamily):
    r"""Gamma distribution using a standard formulation.

    Original formulation presented in [Mourragui et al, 2023].

    family_params of interest:
        - "min_val" (int): min data value, default to 1e-5.
        - "max_val" (int): max data value, default to 10e6.
        - "n_jobs" (int): number of jobs, specifically
        for computing the "nu" parameter.
        - "method" (str): method use to compute the "nu" parameter per
        feature. Two possibles: "MLE" and "MM". Default to "MLE".
        - "eps" (float): minimum difference used for inverting the g
        function. Default to 1e-4
        - "maxiter" (int): maximum number of iterations for the inversion
        of the g function. Default to 100.
    """

    def __init__(self, family_params=None, **kwargs):
        self.family_name = "gamma"
        if family_params is None or "nu" not in family_params:
            print("Gamma distribution not initialized yet")
        default_family_params = {"min_val": 1e-5, "max_val": 10e6, "n_jobs": 1, "eps": 1e-4, "maxiter": 100}
        self.family_params = family_params if family_params else default_family_params
        for k in kwargs:
            self.family_params[k] = kwargs[k]

    def sufficient_statistics(self, X: torch.Tensor):
        return torch.stack([torch.log(X), X])

    def natural_parametrization(self, theta: torch.Tensor):
        nat_params = torch.stack([theta, -torch.ones(theta.shape) * self.family_params["nu"]])
        if nat_params.shape[1] == 1:
            nat_params = nat_params.flatten()
        return nat_params

    def log_partition(self, theta: torch.Tensor):
        return torch.lgamma(theta + 1) - (theta + 1) * torch.log(self.family_params["nu"])

    def exponential_term(self, X: torch.Tensor, theta: torch.Tensor):
        return torch.sum(torch.multiply(self.sufficient_statistics(X), self.natural_parametrization(theta)), axis=0)

    def _digamma_implicit_function(self, X: torch.Tensor, theta: torch.Tensor):
        return torch.digamma(theta + 1) - torch.log(X * self.family_params["nu"])

    def invert_g(self, X: torch.Tensor = None):
        """Dichotomy to compute inverse of digamma function."""

        # Initialize dichotomy parameters.
        min_val = torch.zeros(X.shape)
        max_val = torch.ones(X.shape) * self.family_params["max_val"]
        theta = (min_val + max_val) / 2

        llik = self._digamma_implicit_function(X, theta)
        for idx in tqdm(range(self.family_params["maxiter"])):
            max_val[llik > 0] = theta[llik > 0]
            min_val[llik < 0] = theta[llik < 0]
            theta = (min_val + max_val) / 2
            llik = self._digamma_implicit_function(X, theta)

            if torch.max(torch.abs(llik)) < self.family_params["eps"]:
                print("CONVERGENCE AFTER %s ITERATIONS" % (idx))
                break

        if idx == self.family_params["maxiter"]:
            print("CONVERGENCE NOT REACHED")

        return theta

    def initialize_family_parameters(self, X: torch.Tensor = None):
        p = X.shape[1]

        self.family_params["nu"] = torch.Tensor(
            Parallel(n_jobs=self.family_params["n_jobs"])(
                delayed(scipy.stats.gamma.fit)(X[:, idx], floc=0) for idx in tqdm(range(p))
            )
        )
        self.family_params["nu"] = 1.0 / self.family_params["nu"][:, -1]
        assert self.family_params["nu"].shape[0] == p

        return True


class LogNormal(ExponentialFamily):
    r"""Log-normal distribution using a standard formulation.

    Original formulation presented in [Mourragui et al, 2023].

    family_params of interest:
        - "min_val" (int): min data value, default to 1e-5.
        - "max_val" (int): max data value, default to 10e6.
        - "n_jobs" (int): number of jobs, specifically
        for computing the "nu" parameter.
        - "method" (str): method use to compute the "nu" parameter per
        feature. Two possibles: "MLE" and "MM". Default to "MLE".
        - "eps" (float): minimum difference used for inverting the g
        function. Default to 1e-4
        - "maxiter" (int): maximum number of iterations for the inversion
        of the g function. Default to 100.
    """

    def __init__(self, family_params=None, **kwargs):
        self.family_name = "log_normal"
        if family_params is None or "nu" not in family_params:
            print("Log Normal distribution not initialized yet")
        default_family_params = {
            "min_val": 1e-8,
            "max_val": 10e6,
            "n_jobs": 1,
            "eps": 1e-4,
            "maxiter": 100,
        }
        self.family_params = family_params if family_params else default_family_params
        for k in kwargs:
            self.family_params[k] = kwargs[k]

    def sufficient_statistics(self, X: torch.Tensor):
        log_X = torch.log(X)
        return torch.stack([log_X, torch.square(log_X)])

    def natural_parametrization(self, theta: torch.Tensor):
        nat_params = torch.stack(
            [
                theta / torch.square(self.family_params["nu"]),
                -torch.ones(theta.shape) / (2 * torch.square(self.family_params["nu"])),
            ]
        )
        if nat_params.shape[1] == 1:
            nat_params = nat_params.flatten()
        return nat_params

    def base_measure(self, X: torch.Tensor):
        return 1 / (np.sqrt(2 * torch.pi) * X)

    def log_partition(self, theta: torch.Tensor):
        return torch.square(theta) / (2 * torch.square(self.family_params["nu"])) + torch.log(self.family_params["nu"])

    def exponential_term(self, X: torch.Tensor, theta: torch.Tensor):
        return torch.sum(torch.multiply(self.sufficient_statistics(X), self.natural_parametrization(theta)), axis=0)

    def invert_g(self, X: torch.Tensor = None):
        return torch.log(X.clip(self.family_params["min_val"]))

    def initialize_family_parameters(self, X: torch.Tensor = None):
        self.family_params["nu"] = torch.sqrt(torch.var(torch.log(X), axis=0))
