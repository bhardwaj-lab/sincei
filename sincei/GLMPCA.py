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

from .ExponentialFamily import Normal
from .ExponentialFamily import Poisson
from .ExponentialFamily import Bernoulli

# from .negative_binomial_routines import compute_dispersion
# from .exponential_family import *
# from .log_normal import LOG_NORMAL_ZERO_THRESHOLD


class GLMPCA:
    def __init__(
        self,
        n_pc,
        family,
        maxiter=1000,
        max_param=10,
        learning_rate=0.02,
        batch_size=128,
        step_size=20,
        gamma=0.5,
        n_init=1,
        n_jobs=1,
    ):
        self.n_pc = n_pc
        self.family = family
        self.maxiter = maxiter
        self.log_part_theta_matrices_ = None
        self.max_param = np.abs(max_param)
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
