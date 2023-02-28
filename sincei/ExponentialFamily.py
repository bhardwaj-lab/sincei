import torch
from copy import deepcopy
import numpy as np
import scipy
from joblib import Parallel, delayed

#from .beta_routines import compute_alpha, compute_alpha_gene, compute_mu_gene
#from.log_normal import LOG_NORMAL_ZERO_THRESHOLD, pi_val
#from .gamma import GAMMA_ZERO_THRESHOLD, compute_gamma_saturated_params_gene

saturation_eps = 10**-10

class ExponentialFamily:

    def __init__(self, family_params=None):
        self.family_name = 'base'
        self.family_params = family_params

    def sufficient_statistics(self, X: torch.Tensor):
        return X

    def natural_parametrization(self, theta: torch.Tensor):
        return theta

    def log_partition(self, theta: torch.Tensor):
        return 0.

    def base_measure(self, X: torch.Tensor):
        return 1.

    def invert_g(self, X: torch.Tensor):
        return X

    def exponential_term(self, X: torch.Tensor, theta: torch.Tensor):
        return torch.multiply(
            self.sufficient_statistics(X),
            self.natural_parametrization(theta)
        )

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
        return - torch.sum(
            self.log_distribution(X, theta)
        )

    def compute_ancillary_params(self, X: torch.Tensor=None):
        pass


class Gaussian(ExponentialFamily):
    def __init__(self, family_params=None):
        self.family_name = 'gaussian'
        self.family_params = family_params if family_params else {}

    def sufficient_statistics(self, X: torch.Tensor):
        return X

    def natural_parametrization(self, theta: torch.Tensor):
        return theta

    def log_partition(self, theta: torch.Tensor):
        return torch.square(theta) / 2.

    def base_measure(self, X: torch.Tensor):
        return torch.exp(-torch.square(X)/2.) / np.sqrt(2.*torch.pi)

    def invert_g(self, X: torch.Tensor):
        return X

    def load_family_params_to_gpu(self, device):
        self.family_params = {
            k: self.family_params[k].to(device)
            for j in self.family_params
        }


class Bernoulli(ExponentialFamily):
    def __init__(self, family_params=None):
        self.family_name = 'bernoulli'
        self.family_params = family_params if family_params else {'max_val': 30}

    def sufficient_statistics(self, X: torch.Tensor):
        return X

    def natural_parametrization(self, theta: torch.Tensor):
        return theta

    def log_partition(self, theta: torch.Tensor):
        return torch.log(1.+torch.exp(theta))

    def base_measure(self, X: torch.Tensor):
        return 1.

    def invert_g(self, X: torch.Tensor):
        return torch.log(X/(X-1)).clip(-self.family_params['max_val'], self.family_params['max_val'])


class Poisson(ExponentialFamily):
    def __init(self, family_params=None):
        self.family_name = 'poisson'
        self.family_params = family_params

    def sufficient_statistics(self, X: torch.Tensor):
        return X

    def natural_parametrization(self, theta: torch.Tensor):
        return theta

    def log_partition(self, theta: torch.Tensor):
        return torch.exp(theta)

    def base_measure(self, X: torch.Tensor):
        return 1/scipy.special.gamma(X+1)

    def invert_g(self, X: torch.Tensor):
        return torch.log(X)

    def log_distribution(self, X: torch.Tensor, theta: torch.Tensor):
        """The computation of gamma function for the base measure (h) would lead to inf, hence a re-design
        of the method."""
        log_f = torch.lgamma(X+1)
        expt = self.exponential_term(X, theta) - self.log_partition(theta)

        return expt - log_f


