import numpy as np
import pandas as pd
import torch, os
import torch.optim
from copy import deepcopy
from joblib import Parallel, delayed
import mctorch.nn as mnn
from pickle import dump, load
import mctorch.optim as moptim
from torch.utils.data import Dataset, TensorDataset, DataLoader
from scipy.stats import beta as beta_dst
from scipy.stats import lognorm
from scipy.stats import gamma as gamma_dst

from .ExponentialFamily import Gaussian
from .ExponentialFamily import Poisson
from .ExponentialFamily import Bernoulli

#from .negative_binomial_routines import compute_dispersion
#from .exponential_family import *
#from .log_normal import LOG_NORMAL_ZERO_THRESHOLD

EXPONENTIAL_FAMILY_DICT = {
    'gaussian': Gaussian,
    'poisson': Poisson,
    'binomial': Bernoulli
}

class GLMPCA:

    def __init__(self, n_pc, family,
                 family_params=None,
                 max_iter=1000,
                 learning_rate = 0.02,
                 batch_size=128,
                 step_size=20,
                 gamma=0.5,
                 n_init=1,
                 n_jobs=1):

        self.n_pc = n_pc
        self.family = family
        self.family_params = family_params
        self.log_part_theta_matrices_ = None
        self.max_iter = np.abs(max_iter)
        self.learning_rate_ = learning_rate
        self.initial_learning_rate_ = learning_rate
        self.n_jobs = n_jobs
        self.batch_size = batch_size
        self.n_init = n_init
        self.gamma = gamma
        self.step_size = step_size

        self.saturated_loadings_ = None
        # saturated_intercept_: before projecting
        self.saturated_intercept_ = None
        # reconstruction_intercept: after projecting
        self.reconstruction_intercept_ = None

        # Whether to perform sample or gene projection
        self.sample_projection = False

        self.exp_family_params = None
        self.loadings_learning_scores_ = []
        self.loadings_learning_rates_ = []

        # Set up exponential family
        if type(family) is str:
            self.exponential_family = EXPONENTIAL_FAMILY_DICT[family](self.family_params)
        else:
            self.exponential_family = family

    def fit(self, X):
        # Fit exponential family params (e.g., dispersion for negative binomial)
        self.exponential_family.compute_ancillary_params(X)

        # Compute saturated parameters, alongside exponential family parameters (if needed)
        saturated_parameters = self.exponential_family.invert_g(X)

        # Use saturated parameters to find loadings by projected gradient descent
        self.compute_saturated_loadings(X)

    def compute_saturated_loadings(self, X):
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        # Compute saturated parameters and load on divide
        saturated_parameters = self.exponential_family.invert_g(X)
        saturated_parameters = saturated_parameters.to(device)

        # Initialize the learning procedure
        self.learning_rate_ = self.initial_learning_rate_
        self.loadings_learning_scores_ = []
        self.loadings_learning_rates_ = []

        # Train GLM-PCA with mcTorch.
        loadings_ = self._saturated_loading_iter(saturated_parameters, batch_size=self.batch_size)
        self.saturated_loadings_ = loadings_[0]
        self.saturated_intercept_ = loadings_[1]

    def _saturated_loading_iter(self, saturated_parameters: torch.Tensor, batch_size=128):
        return (None, None)

