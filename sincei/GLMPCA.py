import numpy as np
import torch
import torch.optim
from torch.utils.data import TensorDataset, DataLoader
from tqdm import tqdm
import anndata as ad

try:
    import mctorch.nn as mnn
    import mctorch.optim as moptim
except ImportError:
    raise ImportError("Please install mctorch package via `pip install --user mctorch")

from sincei.ExponentialFamily import Gaussian, Poisson, Bernoulli, Beta, Gamma, LogNormal, SigmoidBeta

EXPONENTIAL_FAMILY_DICT = {
    "gaussian": Gaussian,
    "poisson": Poisson,
    "bernoulli": Bernoulli,
    "beta": Beta,
    "gamma": Gamma,
    "lognormal": LogNormal,
    "log_normal": LogNormal,
    "sigmoid_beta": SigmoidBeta,
}
LEARNING_RATE_LIMIT = 10 ** (-10)


class GLMPCA:
    r"""Performs GLM-PCA on a data matrix to reduce dimensionality

    This class computes the generalized-linear model principal components (GLM-PC)
    of a dataset by exploiting the framework of saturated parameters.
    Specifically, given an exponential family distribution choosen following some
    prior knowledge, GLM-PCA will find a collection of directions which maximise the
    reconstruction error, computed as the negative log-likelihood of the chosen
    exponential family.

    By making use of an alternative formulation, our implementation can exploit
    automatic differentiation and can therefore rely on mini-batch Stochastic
    Gradient Descent. As a consequence, it scales to very large dataset.

    Another interesting feature of our implementation is that it does not require
    cumbersome Lagrangian derivations. If you wish to test an exponential family
    distribution not present in our implementation, adding a class in ExponentialFamily
    with the different density functionals defined there would suffice to use GLM-PCA.

    Parameters
    ----------
    n_pc : int
        Number of principal components.

    family: str
        Name of the exponential family distribution. Possible families: "gaussian", "poisson",
        "bernoulli", "beta", "gamma", "log_normal", "log_beta":, "sigmoid_beta". Default to
        "gaussian"

    family_params : dict
        Dictionary with family parameters to be added. List of parameters depend on the
        specific ExponentialFamily class chosen. Examples:
            - "n_jobs" (int) for parallelization, specifically for "beta" and "gamma".
            - "min_val" (float) for truncating in "poisson" or "beta".
            - "eps" (float) for convergence in inverse computation in "beta".
        Default to None.

    max_iter : int
        Maximum number of epochs in the GLM-PCA optimisation. Default to 100.

    learning_rate: float
        Learning rate to be used in the GLM-PCA optimisation. If learning_rate is too
        high and lead to NaN, our implementation automatically restarts the optimisation
        with a smaller value. Default to 0.2.


    batch_size : int
        Size of the batch in the SGD optimisation step. Default to 256.

    step_size: int
        Step size in optimiser scheduler. See more: https://pytorch.org/docs/stable/generated/torch.optim.lr_scheduler.StepLR.html
        Default to 20.

    gamma: int
        Reduction parameter for optimiser scheduler. See more: https://pytorch.org/docs/stable/generated/torch.optim.lr_scheduler.StepLR.html
        Default to 0.5

    n_init: int
        Number of GLM-PCA initializations. Useful if you want to explore different
        random seeds and starting points. Default to 1.

    init: str
        Method to initialize loadings. "spectral" performs SVD on  the saturated parameters from
        a small random batch of the dataset, "random" performs a random initialization on the
        Stiefel manifold. Default to "spectral".

    n_jobs: int
        Number of jobs used in parallel operations. Default to 1.
    """

    def __init__(
        self,
        n_pc,
        family="gaussian",
        family_params=None,
        max_iter=100,
        learning_rate=0.2,
        batch_size=256,
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
        r"""Fits a GLM-PCA to a specific dataset.

        Parameters
        ----------
        X : torch.Tensor, np.array or anndata object
            Dataset with cells in the rows and features in the columns.

        Returns
        -------
        bool
            Returns True if the fitting procedure has been successful.
        """
        if isinstance(X, ad.AnnData):
            X_fit = torch.Tensor(X.X.transpose().toarray())
        elif isinstance(X, np.ndarray):
            X_fit = torch.Tensor(X)
        elif isinstance(X, torch.Tensor):
            X_fit = X.clone()
        else:
            raise ValueError("X format unrecognised: %s != np.ndarray or torch.Tensor" % (type(X)))

        # Fit exponential family params (e.g., dispersion for negative binomial)
        self.exponential_family.initialize_family_parameters(X_fit)

        # Compute saturated parameters, alongside exponential family parameters (if needed)
        saturated_parameters = self.exponential_family.invert_g(X_fit)

        # Use saturated parameters to find loadings by projected gradient descent
        self.saturated_loadings_ = []
        self.saturated_intercept_ = []

        # Initialize the learning procedure
        self.learning_rate_ = self.initial_learning_rate_
        self.loadings_learning_scores_ = []
        self.loadings_learning_rates_ = []

        # Compute loadings for different parameters
        for _ in range(self.n_init):
            self._compute_saturated_loadings(X_fit, saturated_parameters)

        # Select best model
        training_cost = torch.Tensor(
            [
                self._optim_cost(loadings, intercept, X_fit, saturated_parameters)
                for loadings, intercept in zip(self.saturated_loadings_, self.saturated_intercept_)
            ]
        )
        best_model_idx = torch.argmin(training_cost)
        self.saturated_intercept_ = self.saturated_intercept_[best_model_idx]
        self.saturated_loadings_ = self.saturated_loadings_[best_model_idx]

        return True

    def transform(self, X):
        r"""Transforms and project dataset X onto the principal components.

        Parameters
        ----------
        X : torch.Tensor or np.array
            Dataset with cells in the rows and features in the columns.

        Returns
        -------
        torch.Tensor
            Projected saturated parameters.
        """
        saturated_parameters = self.exponential_family.invert_g(X)

        # Compute intercept term
        n = X.shape[0]
        intercept_term = self.saturated_intercept_.unsqueeze(0).repeat(n, 1).to(self.device)

        projected_parameters = saturated_parameters - intercept_term
        projected_parameters = projected_parameters.matmul(self.saturated_loadings_)

        return projected_parameters

    def _compute_saturated_loadings(self, X, saturated_parameters):
        # Runs one optimisation of the loadings, using mcTorch.
        loadings_ = self._saturated_loading_iter(saturated_parameters, X)
        self.saturated_loadings_.append(loadings_[0])
        self.saturated_intercept_.append(loadings_[1])

    def _saturated_loading_iter(self, saturated_parameters: torch.Tensor, X: torch.Tensor):
        r"""Computes the loadings solution of the GLM-PCA optimisation problem.

        Parameters
        ----------
        saturated_parameters : torch.Tensor
            Saturated parameters of the dataset X ($g^{-1}\left(X\right)$)

        X : torch.Tensor
            Dataset with cells in the rows and features in the columns.

        Returns
        -------
        torch.Tensor
            Projected saturated parameters.
        """
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

        # Load dataset
        train_data = TensorDataset(X, saturated_parameters.data.clone())
        train_loader = DataLoader(dataset=train_data, batch_size=self.batch_size, shuffle=True, drop_last=True)

        # Run epoch in a for loop
        self._loadings_epochs = [_loadings.clone().detach()]
        self._intercept_epochs = [_intercept.clone().detach()]
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

            self._loadings_epochs.append(_loadings.clone().detach())
            self._intercept_epochs.append(_intercept.clone().detach())

            # If NaN or Inf is found in the parameters, start over optimisation with reduced learning rate.
            if np.isinf(self.loadings_learning_scores_[-1][-1]) or np.isnan(self.loadings_learning_scores_[-1][-1]):
                print("\tRESTART BECAUSE INF/NAN FOUND", flush=True)
                self.learning_rate_ = self.learning_rate_ * self.gamma
                self.loadings_learning_scores_ = self.loadings_learning_scores_[:-1]
                self.loadings_learning_rates_ = self.loadings_learning_rates_[:-1]

                # Remove memory
                del train_data, train_loader, _optimizer, _loadings, _intercept, _lr_scheduler
                if "cuda" in str(self.device):
                    torch.cuda.empty_cache()

                self._loadings_epochs = []
                self._intercept_epochs = []

                return self._saturated_loading_iter(
                    saturated_parameters=saturated_parameters,
                    X=X,
                )

        return (_loadings, _intercept)

    def _create_saturated_loading_optim(self, parameters: torch.Tensor, X: torch.Tensor):
        r"""Initialises the optimisation problem.

        Parameters
        ----------
        saturated_parameters : torch.Tensor
            Saturated parameters of the dataset X ($g^{-1}\left(X\right)$)

        X : torch.Tensor
            Dataset with cells in the rows and features in the columns.

        Returns
        -------
        optimizer: mctorch.optim
            mctorch optimiser instance

        loadings: mnn.Parameter
            mcTorch parameter instance with loadings constrained to the Stiefel manifold

        intercept: mnn.Parameter
            mcTorch parameter instance with intercept.

        lr_scheduler: torch.optim.scheduler
            Scheduler instance.
        """
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        # Initialize loadings with spectrum (2**13 as maximum value for SVD to be relatively fast)
        random_batch_size = min(X.shape[0], 2**13)
        random_idx = np.random.choice(np.arange(parameters.shape[0]), replace=False, size=random_batch_size)
        if self.init == "spectral":
            _, _, v = torch.linalg.svd(parameters[random_idx] - torch.mean(parameters[random_idx], axis=0))
            loadings = mnn.Parameter(
                data=v[: self.n_pc, :].T, manifold=mnn.Stiefel(parameters.shape[1], self.n_pc), requires_grad=True
            )  # .to(
            # self.device
            # )
        elif self.init == "random":
            loadings = mnn.Parameter(
                manifold=mnn.Stiefel(parameters.shape[1], self.n_pc), requires_grad=True
            )  # .to(self.device)

        # Initialize intercept
        if self.exponential_family.family_name in ["poisson"]:
            intercept = mnn.Parameter(
                torch.median(parameters[random_idx], axis=0).values,
                manifold=mnn.Euclidean(parameters.shape[1]),
                requires_grad=True,
            )
            # ).to(self.device)
        else:
            intercept = mnn.Parameter(
                torch.mean(parameters[random_idx], axis=0),
                manifold=mnn.Euclidean(parameters.shape[1]),
                requires_grad=True,
            )
            # ).to(self.device)

        # Load ExponentialFamily params to GPU (if they exist)
        self.exponential_family.load_family_params_to_gpu(self.device)

        # Create optimizer
        # TODO: allow for other optimizer to be used.
        # TODO: learning rate for intercept.
        print("LEARNING RATE: %s" % (self.learning_rate_))
        optimizer = moptim.rAdagrad(
            params=[
                {"params": loadings, "lr": self.learning_rate_},
                {"params": intercept, "lr": self.learning_rate_ * 0.01},
            ]
        )
        lr_scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=self.step_size, gamma=self.gamma)

        return optimizer, loadings, intercept, lr_scheduler

    def _optim_cost(
        self, loadings: torch.Tensor, intercept: torch.Tensor, batch_data: torch.Tensor, batch_parameters: torch.Tensor
    ):
        n = batch_data.shape[0]
        intercept_term = intercept.unsqueeze(0).repeat(n, 1)  # .to(self.device)

        projected_parameters = batch_parameters - intercept_term
        projected_parameters = projected_parameters.matmul(loadings).matmul(loadings.T)
        projected_parameters = projected_parameters + intercept_term

        return self.exponential_family.log_likelihood(batch_data, projected_parameters)
