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
from tqdm import tqdm
from scipy.stats import beta as beta_dst
from scipy.stats import lognorm
from scipy.stats import gamma as gamma_dst

from sincei.ExponentialFamily import Gaussian, Poisson, Bernoulli, Beta, Gamma, LogNormal

EXPONENTIAL_FAMILY_DICT = {
    "gaussian": Gaussian,
    "poisson": Poisson,
    "bernoulli": Bernoulli,
    "beta": Beta,
    "gamma": Gamma,
    "lognormal": LogNormal,
    "log_normal": LogNormal,
}
LEARNING_RATE_LIMIT = 10 ** (-10)


class GLMPCA:
    def __init__(
        self,
        n_pc,
        family,
        family_params=None,
        max_iter=1000,
        learning_rate=0.02,
        batch_size=128,
        step_size=20,
        gamma=0.5,
        n_init=1,
        init="spectral",
        n_jobs=1,
    ):
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
        self.init = init

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

        # Initialize device
        self.device = None

        # Set up exponential family
        if type(family) is str:
            self.exponential_family = EXPONENTIAL_FAMILY_DICT[family](self.family_params)
        else:
            self.exponential_family = family

    def fit(self, X):
        # Fit exponential family params (e.g., dispersion for negative binomial)
        self.exponential_family.initialize_family_parameters(X)

        # Compute saturated parameters, alongside exponential family parameters (if needed)
        saturated_parameters = self.exponential_family.invert_g(X)

        # Use saturated parameters to find loadings by projected gradient descent
        self.saturated_loadings_ = []
        self.saturated_intercept_ = []

        # Initialize the learning procedure
        self.learning_rate_ = self.initial_learning_rate_
        self.loadings_learning_scores_ = []
        self.loadings_learning_rates_ = []

        # Compute loadings for different parameters
        for _ in range(self.n_init):
            self._compute_saturated_loadings(X)

        # Select best model
        training_cost = torch.Tensor(
            [
                self._optim_cost(loadings, intercept, X, saturated_parameters)
                for loadings, intercept in zip(self.saturated_loadings_, self.saturated_intercept_)
            ]
        )
        best_model_idx = torch.argmin(training_cost)
        self.saturated_intercept_ = self.saturated_intercept_[best_model_idx]
        self.saturated_loadings_ = self.saturated_loadings_[best_model_idx]

    def transform(self, X):
        saturated_parameters = self.exponential_family.invert_g(X)

        # Compute intercept term
        n = X.shape[0]
        intercept_term = self.saturated_intercept_.unsqueeze(0).repeat(n, 1).to(self.device)

        projected_parameters = saturated_parameters - intercept_term
        projected_parameters = projected_parameters.matmul(self.saturated_loadings_)

        return projected_parameters

    def _compute_saturated_loadings(self, X):
        # Compute saturated parameters and load on divide
        saturated_parameters = self.exponential_family.invert_g(X)

        # Train GLM-PCA with mcTorch.
        loadings_ = self._saturated_loading_iter(saturated_parameters, X)
        self.saturated_loadings_.append(loadings_[0])
        self.saturated_intercept_.append(loadings_[1])

    def _saturated_loading_iter(self, saturated_parameters: torch.Tensor, X: torch.Tensor):
        if self.learning_rate_ < LEARNING_RATE_LIMIT:
            raise ValueError("LEARNING RATE IS TOO SMALL : DID NOT CONVERGE")

        # Set device for GPU usage
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        # Set up list of saving scores
        self.loadings_learning_scores_.append([])
        self.loadings_learning_rates_.append([])

        _optimizer, _loadings, _intercept, _lr_scheduler = self._create_saturated_loading_optim(
            parameters=saturated_parameters.data.clone(), X=X
        )

        train_data = TensorDataset(X, saturated_parameters.data.clone())
        train_loader = DataLoader(dataset=train_data, batch_size=self.batch_size, shuffle=True)

        for _ in tqdm(range(self.max_iter)):
            loss_val = []
            for batch_data, batch_parameters in train_loader:
                cost_step = self._optim_cost(
                    loadings=_loadings,
                    intercept=_intercept,
                    batch_data=batch_data,
                    batch_parameters=batch_parameters,
                )

                if "cuda" in str(self.device):
                    self.loadings_learning_scores_[-1].append(cost_step.cpu().detach().numpy())
                else:
                    self.loadings_learning_scores_[-1].append(cost_step.detach().numpy())
                cost_step.backward()
                _optimizer.step()
                _optimizer.zero_grad()
                self.loadings_learning_rates_[-1].append(_lr_scheduler.get_last_lr())
            _lr_scheduler.step()

            if np.isinf(self.loadings_learning_scores_[-1][-1]) or np.isnan(self.loadings_learning_scores_[-1][-1]):
                print("\tRESTART BECAUSE INF/NAN FOUND", flush=True)
                self.learning_rate_ = self.learning_rate_ * self.gamma
                self.loadings_learning_scores_ = self.loadings_learning_scores_[:-1]
                self.loadings_learning_rates_ = self.loadings_learning_rates_[:-1]

                # Remove memory
                del train_data, train_loader, _optimizer, _loadings, _intercept, _lr_scheduler
                if "cuda" in str(self.device):
                    torch.cuda.empty_cache()

                return self._saturated_loading_iter(
                    saturated_parameters=saturated_parameters,
                    X=X,
                )

        return (_loadings, _intercept)

    def _create_saturated_loading_optim(self, parameters: torch.Tensor, X: torch.Tensor):
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        # Initialize loadings with spectrum
        random_idx = np.random.choice(np.arange(parameters.shape[0]), replace=False, size=self.batch_size)
        if self.init == "spectral":
            _, _, v = torch.linalg.svd(parameters[random_idx] - torch.mean(parameters[random_idx], axis=0))
            loadings = mnn.Parameter(data=v[: self.n_pc, :].T, manifold=mnn.Stiefel(parameters.shape[1], self.n_pc)).to(
                self.device
            )
        elif self.init == "random":
            loadings = mnn.Parameter(manifold=mnn.Stiefel(parameters.shape[1], self.n_pc)).to(self.device)

        # Initialize intercept
        if self.exponential_family.family_name in ["poisson"]:
            intercept = mnn.Parameter(
                torch.median(parameters[random_idx], axis=0).values, manifold=mnn.Euclidean(parameters.shape[1])
            ).to(self.device)
        else:
            intercept = mnn.Parameter(
                torch.mean(parameters[random_idx], axis=0), manifold=mnn.Euclidean(parameters.shape[1])
            ).to(self.device)

        # Load ExponentialFamily params to GPU (if they exist)
        self.exponential_family.load_family_params_to_gpu(self.device)

        # Create optimizer
        # TODO: allow for other optimizer to be used.
        print("LEARNING RATE: %s" % (self.learning_rate_))
        optimizer = moptim.rAdagrad(
            params=[{"params": loadings, "lr": self.learning_rate_}, {"params": intercept, "lr": self.learning_rate_}]
        )
        lr_scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=self.step_size, gamma=self.gamma)

        return optimizer, loadings, intercept, lr_scheduler

    def _optim_cost(
        self, loadings: torch.Tensor, intercept: torch.Tensor, batch_data: torch.Tensor, batch_parameters: torch.Tensor
    ):
        n = batch_data.shape[0]
        intercept_term = intercept.unsqueeze(0).repeat(n, 1).to(self.device)

        projected_parameters = batch_parameters - intercept_term
        projected_parameters = projected_parameters.matmul(loadings).matmul(loadings.T)
        projected_parameters = projected_parameters + intercept_term

        return self.exponential_family.log_likelihood(batch_data, projected_parameters)
