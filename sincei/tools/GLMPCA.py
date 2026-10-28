from __future__ import annotations

import copy
import time
from dataclasses import dataclass
from typing import Any, Literal

import anndata as ad
import numpy as np
import torch
from scipy import sparse as sp
from tqdm import tqdm

from sincei.tools.ExponentialFamily import ExponentialFamily, GLMFamily

try:
    from anndata.abc import CSRDataset
except ImportError:  # anndata < 0.11
    try:
        from anndata.experimental import CSRDataset
    except ImportError:
        CSRDataset = ()  # type: ignore[assignment]

LEARNING_RATE_LIMIT = 1e-10
# a densified chunk holds at most this many values (~200 MB in float32)
MAX_CHUNK_ENTRIES = 5 * 10**7
# families supporting total-count depth correction (normalize: X / s; offset: log s on theta)
DEPTH_CORRECTION_FAMILIES = {"normalize": ("poisson", "negative_binomial", "gamma", "lognormal"),
                             "offset": ("poisson", "negative_binomial", "lognormal")}


# =====================================================================================
# Data access and saturated parameters
# =====================================================================================
@dataclass
class _Batch:
    X: torch.Tensor | None = None      # dense data (optionally a column subset)
    theta: torch.Tensor | None = None  # dense saturated parameters (dense input)
    delta: torch.Tensor | None = None  # sparse Theta - 1 theta0^T (sparse input)
    o: torch.Tensor | None = None      # per-cell log size factor (depth_correction="offset")


class _Data:
    r"""Chunked row access and saturated parameters for a set of cells.

    Accepts numpy / torch / scipy.sparse / AnnData (also backed='r'); the full
    matrix is never densified. For sparse X, g^{-1} is element-wise, hence
    Theta = 1 theta0^T + Delta with theta0 = g^{-1}(0) per feature and Delta sparse on
    the non-zeros of X: scores cost O(nnz d) and the cache O(nnz). Exact for any
    family.
    """

    def __init__(self, X, device, dtype, chunk_size: int, transpose: bool = False,
                 n_sub: int | None = None, rng: np.random.Generator | None = None,
                 size_factors: np.ndarray | None = None, depth_correction: str | None = None):
        self.s, self.depth_correction = size_factors, depth_correction  # s indexed by global row
        backed = isinstance(X, ad.AnnData) and X.isbacked
        M = X.X if isinstance(X, ad.AnnData) else X
        if isinstance(M, torch.Tensor):
            M = M.detach().cpu().numpy()
        if transpose:
            if backed:
                raise ValueError("sample_projection=True requires in-memory data.")
            M = M.T
        if sp.issparse(M):
            M = sp.csr_matrix(M)
            M.sum_duplicates()
        self.M, self.device, self.dtype = M, device, dtype
        self.sparse = sp.issparse(M) or isinstance(M, CSRDataset)
        n_all, self.p = M.shape
        self.chunk = max(1, min(chunk_size, MAX_CHUNK_ENTRIES // self.p))
        self.cells = np.arange(n_all)
        if n_sub is not None and n_sub < n_all:
            self.cells = np.sort(rng.choice(n_all, n_sub, replace=False))
        self.n = len(self.cells)
        # rows held in memory: in-memory input, or a subsample read once from disk
        if not backed:
            self.mem = M if self.n == n_all else M[self.cells]
            if depth_correction == "normalize":
                self.mem = self._scale(self.mem, self.cells)
        else:
            self.mem = self._read(self.cells) if self.n < n_all else None
        self.family, self.theta0, self.cache = None, None, None

    # -- rows ------------------------------------------------------------------------
    def _read(self, rows: np.ndarray):
        out = self.M[rows]
        out = sp.csr_matrix(out) if self.sparse else np.asarray(out)
        return self._scale(out, rows) if self.depth_correction == "normalize" else out

    def _scale(self, Xc, rows: np.ndarray):
        """X / s row-wise (CSR stays CSR)."""
        inv = (1.0 / self.s[rows]).astype(np.float32)
        return sp.csr_matrix(sp.diags(inv) @ Xc) if sp.issparse(Xc) else np.asarray(Xc) * inv[:, None]

    def rows(self, idx: np.ndarray):
        """Raw rows (CSR or ndarray) for sorted local indices."""
        return self.mem[idx] if self.mem is not None else self._read(self.cells[idx])

    def dense(self, Xc, cols: torch.Tensor | None = None) -> torch.Tensor:
        if cols is not None:
            Xc = Xc[:, cols.cpu().numpy()]
        Xc = Xc.toarray() if sp.issparse(Xc) else np.asarray(Xc)
        return torch.from_numpy(Xc).to(self.device, self.dtype)

    def chunks(self, idx: np.ndarray) -> list[np.ndarray]:
        return [idx[s:s + self.chunk] for s in range(0, len(idx), self.chunk)]

    def minibatches(self, batch_size: int, rng: np.random.Generator) -> list[np.ndarray]:
        if batch_size >= self.n:
            return [np.arange(self.n)]
        if self.mem is None:  # streaming from disk: shuffle contiguous blocks
            starts = rng.permutation(np.arange(0, self.n, batch_size))
            return [np.arange(s, min(s + batch_size, self.n)) for s in starts]
        perm = rng.permutation(self.n)
        return [np.sort(perm[s:s + batch_size]) for s in range(0, self.n, batch_size)]

    # -- saturated parameters --------------------------------------------------------
    def set_family(self, family: ExponentialFamily, cache: bool | str = "auto",
                   max_cache_gb: float = 4.0) -> None:
        self.family = family
        zeros = torch.zeros(1, self.p, device=self.device, dtype=self.dtype)
        self.theta0 = family.invert_g(zeros).reshape(-1)
        if cache is False or self.mem is None:
            return
        if self.sparse:
            # Delta has the sparsity pattern of X: store its values with X's indices
            X, step = self.mem, MAX_CHUNK_ENTRIES
            vals = [self._delta_values(X.data[s:s + step], X.indices[s:s + step]).cpu()
                    for s in range(0, X.nnz, step)]
            vals = torch.cat(vals).numpy() if vals else np.zeros(0)
            self.cache = sp.csr_matrix((vals, X.indices, X.indptr), shape=X.shape)
        elif cache is True or self.n * self.p * torch.finfo(self.dtype).bits / 8 <= max_cache_gb * 1e9:
            self.cache = torch.cat([self.batch(i, need_X=False).theta.cpu()
                                    for i in self.chunks(np.arange(self.n))])

    def _delta_values(self, data: np.ndarray, cols: np.ndarray) -> torch.Tensor:
        cols = torch.from_numpy(cols.astype(np.int64)).to(self.device)
        vals = torch.from_numpy(np.asarray(data)).to(self.device, self.dtype)
        return self.family.subset_features(cols).invert_g(vals) - self.theta0[cols]

    def _to_torch_sparse(self, Xc: sp.csr_matrix, values: torch.Tensor) -> torch.Tensor:
        coo = Xc.tocoo()
        idx = torch.from_numpy(np.vstack([coo.row, coo.col]).astype(np.int64))
        return torch.sparse_coo_tensor(idx.to(self.device), values, Xc.shape, check_invariants=False).coalesce()

    def batch(self, idx: np.ndarray, need_X: bool = True,
              cols: torch.Tensor | None = None) -> _Batch:
        """Data and saturated parameters of the rows ``idx`` (sorted, local)."""
        Xc, b = self.rows(idx), _Batch()
        if self.depth_correction == "offset":
            b.o = torch.as_tensor(np.log(self.s[self.cells[idx]]), device=self.device, dtype=self.dtype)
        if need_X:
            b.X = self.dense(Xc, cols)
        if self.sparse:
            if self.cache is not None:
                D = self.cache[idx]
                vals = torch.from_numpy(D.data).to(self.device, self.dtype)
                b.delta = self._to_torch_sparse(D, vals)
            else:
                b.delta = self._to_torch_sparse(Xc, self._delta_values(Xc.data, Xc.indices))
        elif self.cache is not None:
            b.theta = self.cache[torch.from_numpy(idx)].to(self.device)
        else:
            b.theta = self.family.invert_g(b.X if need_X and cols is None else self.dense(Xc))
        return b

    def scores(self, V: torch.Tensor, mu: torch.Tensor, b: _Batch) -> torch.Tensor:
        """U = (Theta_b - o [x != 0] - 1 mu^T) V: with an offset, non-zeros are shifted by
        -o_i (= g^{-1}(x / s) for log-link families) and zeros stay at theta0."""
        if b.delta is not None:
            U = (self.theta0 - mu) @ V + torch.sparse.mm(b.delta, V)
            nz = torch.sparse_coo_tensor(b.delta.indices(), torch.ones_like(b.delta.values()), b.delta.shape)
        else:
            U, nz = (b.theta - mu) @ V, None
        if b.o is None:
            return U
        nzV = torch.sparse.mm(nz, V) if nz is not None else (b.theta != self.theta0).to(V.dtype) @ V
        return U - b.o[:, None] * nzV


# =====================================================================================
class GLMPCA:
    r"""Performs GLM-PCA on a data matrix to reduce its dimensionality.

    Loadings V (Stiefel) and intercept mu minimise
    L((Theta - 1 mu^T) V V^T + 1 mu^T ; X) with Theta = g^{-1}(X)
    (Landgraf & Lee, 2020; Mourragui et al., 2023). New distributions are added by
    subclassing ``sincei.tools.ExponentialFamily.ExponentialFamily``.

    Parameters
    ----------
    n_pc, family, family_params, max_iter, batch_size, n_init : as before.
    learning_rate : float
        Riemannian Adam step for the loadings (geoopt), scaled internally by
        1/sqrt(n_features) since Stiefel entries are O(1/sqrt(p)).
    gamma : float
        Learning-rate decay factor on plateaus and after a NaN restart.
    init : "spectral", "random" or array (n_features, n_pc) for warm starts.
    optimizer : {"geoopt", "lbfgs"}
        "geoopt": mini-batch Riemannian Adam (GPU, many cells).
        "lbfgs": deterministic full-batch L-BFGS (strong-Wolfe line search) on
        V = qf(W), the QR retraction of an unconstrained W; exact chunked gradients,
        no learning rate. Recommended up to ~1e5 cells, or with fit_mode="subsample".
    fit_mode : {"exact", "minibatch", "subsample"}
        "exact": full objective, chunked (geoopt: one step per epoch).
        "minibatch": stochastic mini-batches of cells (geoopt only).
        "subsample": fit on ``n_cells_fit`` random cells; ``transform`` projects all.
    feature_subsample : float, int or None
        Fraction or number of features at which the likelihood is evaluated per step
        (unbiased). Scores always use all features. Fixed subset with lbfgs.
    n_cells_fit : int or None
        Cells used with fit_mode="subsample" (default min(n, 200_000)).
    chunk_size : int
        Rows per chunk for full passes; capped so a dense chunk <= MAX_CHUNK_ENTRIES.
    intercept_lr_ratio : float
    tol, patience, min_epochs, n_lr_decay :
        Plateau handling (geoopt), via ``torch.optim.lr_scheduler.ReduceLROnPlateau``
        on the epoch NLL per matrix entry: if it improves by less than ``tol`` (abs.)
        for ``patience`` epochs, the learning rate is multiplied by ``gamma``; the
        fit stops at the (n_lr_decay + 1)-th plateau.
    cache_saturated : {"auto", True, False}; max_cache_gb : float
        Keep saturated parameters of in-memory fitting cells (sparse: O(nnz);
        dense: if below max_cache_gb, or always with True).
    whiten : bool
        Rotate V within its span so scores are uncorrelated and ordered (Def. 2.8).
    center_scores : bool
        With whiten, centre the scores via mu <- mu + V mean(U) (likelihood unchanged).
    depth_correction : {None, "normalize", "offset"}
        Per-cell total-count correction, s_i = libsize_i / median(libsize):
        "normalize" fits X / s (poisson, negative_binomial, gamma, lognormal); "offset"
        keeps X in the likelihood and adds log s_i to theta (poisson, negative_binomial,
        lognormal). Size factors can also be passed to fit / transform / score.
    sample_projection : bool
        Transpose the matrix (features become samples).
    device, dtype, seed, verbose, store_epochs.
    """

    def __init__(
        self,
        n_pc: int,
        family: str | ExponentialFamily = "gaussian",
        family_params: dict[str, Any] | None = None,
        max_iter: int = 100,
        learning_rate: float = 0.5,
        batch_size: int = 1024,
        gamma: float = 0.5,
        n_init: int = 1,
        init: Literal["spectral", "random"] | np.ndarray = "spectral",
        optimizer: Literal["geoopt", "lbfgs"] = "geoopt",
        fit_mode: Literal["exact", "minibatch", "subsample"] = "minibatch",
        feature_subsample: float | int | None = None,
        n_cells_fit: int | None = None,
        chunk_size: int = 4096,
        intercept_lr_ratio: float = 1.0,
        tol: float = 1e-5,
        patience: int = 3,
        min_epochs: int = 10,
        n_lr_decay: int = 3,
        cache_saturated: bool | Literal["auto"] = "auto",
        max_cache_gb: float = 4.0,
        whiten: bool = True,
        center_scores: bool = False,
        depth_correction: Literal["normalize", "offset"] | None = None,
        sample_projection: bool = False,
        device: str | torch.device | None = None,
        dtype: torch.dtype = torch.float32,
        seed: int | None = None,
        verbose: bool = True,
        store_epochs: bool = False,
    ) -> None:
        if optimizer not in ("geoopt", "lbfgs"):
            raise ValueError("optimizer must be 'geoopt' or 'lbfgs'")
        if fit_mode not in ("exact", "minibatch", "subsample"):
            raise ValueError("fit_mode must be 'exact', 'minibatch' or 'subsample'")
        if optimizer == "lbfgs" and fit_mode == "minibatch":
            raise ValueError("lbfgs is full-batch: use fit_mode='exact' or "
                             "'subsample', or optimizer='geoopt'.")
        self.n_pc, self.family, self.family_params = n_pc, family, family_params
        self.max_iter, self.batch_size, self.gamma = int(abs(max_iter)), batch_size, gamma
        self.initial_learning_rate_ = self.learning_rate_ = learning_rate
        self.n_init, self.init, self.optimizer, self.fit_mode = n_init, init, optimizer, fit_mode
        self.feature_subsample, self.n_cells_fit = feature_subsample, n_cells_fit
        self.chunk_size = chunk_size
        self.intercept_lr_ratio, self.tol, self.patience = intercept_lr_ratio, tol, patience
        self.min_epochs, self.n_lr_decay = min_epochs, n_lr_decay
        self.cache_saturated, self.max_cache_gb = cache_saturated, max_cache_gb
        self.whiten, self.sample_projection = whiten, sample_projection
        self.center_scores, self.depth_correction = center_scores, depth_correction
        self.dtype, self.seed, self.verbose, self.store_epochs = dtype, seed, verbose, store_epochs
        self.device = torch.device(device or ("cuda" if torch.cuda.is_available() else "cpu"))

        self.saturated_loadings_: torch.Tensor | None = None
        self.saturated_intercept_: torch.Tensor | None = None
        self.reconstruction_intercept_ = self.singular_values_ = self.fit_time_ = None
        self.loadings_learning_scores_: list[list[float]] = []
        self.loadings_learning_rates_: list[list[float]] = []

        if isinstance(family, str):
            family = GLMFamily(family).distribution()(family_params)
        if depth_correction is not None and (
                family.family_name not in DEPTH_CORRECTION_FAMILIES.get(depth_correction, ())
                or sample_projection):
            raise ValueError(f"depth_correction={depth_correction!r} is not supported for "
                             f"{family.family_name!r} (or with sample_projection).")
        # fit() estimates nuisance parameters on a copy, so refits start clean
        self._family_template = self.exponential_family = family

    def _log(self, msg: str) -> None:
        if self.verbose:
            print(msg, flush=True)

    def _feature_sample(self, p: int) -> torch.Tensor | None:
        fs = self.feature_subsample
        if fs is None or self.fit_mode == "exact":
            return None
        k = int(round(fs * p)) if isinstance(fs, float) else int(fs)
        if k >= p:
            return None
        cols = np.sort(self._rng.choice(p, max(k, 1), replace=False))
        return torch.from_numpy(cols).to(self.device)

    def _size_factors(self, X, size_factors=None, train: bool = False) -> np.ndarray | None:
        """s_i = libsize_i / median(training libsize), or the user's size factors."""
        if self.depth_correction is None:
            return None
        if size_factors is None:
            raw = _Data(X, self.device, self.dtype, self.chunk_size)
            lib = np.concatenate([np.asarray(raw.rows(c).sum(1)).ravel()
                                  for c in raw.chunks(np.arange(raw.n))])
            if train:
                self.lib_median_ = float(np.median(lib[lib > 0]))
            size_factors = lib / self.lib_median_
        s = np.asarray(size_factors, dtype=np.float64).ravel()
        return np.maximum(s, s[s > 0].min())  # empty cells: smallest positive factor

    def _data(self, X, n_sub: int | None = None, cache: bool | str = False,
              size_factors=None) -> _Data:
        data = _Data(X, self.device, self.dtype, self.chunk_size, self.sample_projection,
                     n_sub, getattr(self, "_rng", None),
                     self._size_factors(X, size_factors), self.depth_correction)
        data.set_family(self.exponential_family, cache, self.max_cache_gb)
        return data

    # ---------------------------------------------------------------------------------
    def fit(self, X: torch.Tensor | np.ndarray | sp.spmatrix | ad.AnnData,
            size_factors: np.ndarray | None = None) -> bool:
        r"""Fits GLM-PCA on cells (rows) x features (columns). Accepts dense arrays,
        tensors, scipy sparse matrices and AnnData (including backed='r')."""
        t0 = time.time()
        self._rng = np.random.default_rng(self.seed)
        if self.seed is not None:
            torch.manual_seed(self.seed)
        fam = self.exponential_family = copy.deepcopy(self._family_template)
        n_sub = (self.n_cells_fit or 200_000) if self.fit_mode == "subsample" else None
        self.size_factors_ = self._size_factors(X, size_factors, train=True)
        data = _Data(X, self.device, self.dtype, self.chunk_size, self.sample_projection,
                     n_sub, self._rng, self.size_factors_, self.depth_correction)
        if self.n_pc >= data.p:
            raise ValueError("n_pc must be smaller than the number of features.")

        # streamed column statistics: nuisance parameters and column means
        idx = np.arange(data.n)
        if data.n > 500_000:
            idx = np.sort(self._rng.choice(data.n, 500_000, replace=False))
        colsum = torch.zeros(data.p, dtype=torch.float64, device=self.device)
        stats: dict[str, torch.Tensor] = {}
        n_eff = float(len(idx))  # offset: column mean of X / s, i.e. sum x / sum s
        if self.depth_correction == "offset":
            n_eff = float(self.size_factors_[data.cells[idx]].sum())
        for sub in data.chunks(idx):
            Xd = data.dense(data.rows(sub))
            colsum += Xd.double().sum(0)
            if self.depth_correction == "offset":  # nuisance parameters from X / s
                Xd = Xd / torch.as_tensor(self.size_factors_[data.cells[sub]]).to(Xd)[:, None]
            if fam.needs_family_statistics():
                for k, v in fam.family_statistics(Xd).items():
                    stats[k] = stats.get(k, 0) + v
        if fam.needs_family_statistics():
            fam.initialize_family_parameters_from_statistics(stats, len(idx))
        fam.to(self.device, self.dtype)
        self._colmean = (colsum / n_eff).to(self.device, self.dtype)
        data.set_family(fam, self.cache_saturated, self.max_cache_gb)
        self._log(f"GLM-PCA: n={data.n}, p={data.p}, d={self.n_pc}, "
                  f"family={fam.family_name}, optimizer={self.optimizer}, "
                  f"mode={self.fit_mode}, sparse={data.sparse}, "
                  f"cached={data.cache is not None}, device={self.device}")

        self.loadings_learning_scores_, self.loadings_learning_rates_ = [], []
        runs = [self._saturated_loading_iter(data) for _ in range(self.n_init)]

        # select the best run on (up to 50k of) the fitting cells
        idx = np.arange(data.n)
        if data.n > 50_000:
            idx = np.sort(self._rng.choice(data.n, 50_000, replace=False))
        with torch.no_grad():
            costs = [self._accumulate(V, mu, data, idx, scale=1 / len(idx)) for V, mu in runs]
        V, mu = runs[int(np.argmin(costs))]
        self.training_cost_ = min(costs)

        if self.whiten:
            V, self.singular_values_, mu = self._whiten(V, mu, data)
        self.saturated_loadings_, self.saturated_intercept_ = V.detach(), mu.detach()
        self.reconstruction_intercept_ = self.saturated_intercept_
        self.fit_time_ = time.time() - t0
        self._log(f"GLM-PCA fitted in {self.fit_time_:.1f}s, "
                  f"cost/cell={self.training_cost_:.4f}")
        return True

    def _check_fitted(self) -> tuple[torch.Tensor, torch.Tensor]:
        if self.saturated_loadings_ is None:
            raise RuntimeError("GLMPCA is not fitted. Call fit() first.")
        return self.saturated_loadings_, self.saturated_intercept_

    @torch.no_grad()
    def transform(self, X, chunk_size: int | None = None, size_factors=None) -> torch.Tensor:
        r"""Projects X onto the components, (g^{-1}(X) - 1 mu^T) V, chunk by chunk
        (exact out-of-sample extension). Returns a (n_cells, n_pc) CPU tensor."""
        V, mu = self._check_fitted()
        data = self._data(X, size_factors=size_factors)
        if chunk_size:
            data.chunk = chunk_size
        if data.p != V.shape[0]:
            raise ValueError(f"X has {data.p} features, model has {V.shape[0]}.")
        return torch.cat([data.scores(V, mu, data.batch(i, need_X=False)).cpu()
                          for i in data.chunks(np.arange(data.n))])

    @torch.no_grad()
    def score(self, X, size_factors=None) -> float:
        r"""Mean negative log-likelihood per cell of X under the out-of-sample
        reconstruction: a held-out criterion for n_pc, truncation values, etc."""
        V, mu = self._check_fitted()
        data = self._data(X, size_factors=size_factors)
        return self._accumulate(V, mu, data, np.arange(data.n), scale=1 / data.n)

    def fit_transform(self, X, size_factors=None) -> torch.Tensor:
        self.fit(X, size_factors)
        return self.transform(X, size_factors=size_factors)

    # ---------------------------------------------------------------------------------
    # Cost
    # ---------------------------------------------------------------------------------
    def _nll(self, V, mu, data: _Data, b: _Batch, cols=None) -> torch.Tensor:
        """Batch NLL; with ``cols`` evaluated on a feature subset, rescaled by p/|cols|."""
        U = data.scores(V, mu, b)
        o = 0.0 if b.o is None else b.o[:, None]
        if cols is None:
            return self.exponential_family.log_likelihood(b.X, o + mu + U @ V.T)
        fam = self.exponential_family.subset_features(cols)
        return fam.log_likelihood(b.X, o + mu[cols] + U @ V[cols].T) * (data.p / len(cols))

    def _accumulate(self, V, mu, data: _Data, idx: np.ndarray, cols=None,
                    backward: bool = False, scale: float = 1.0) -> float:
        """scale * NLL over rows ``idx``, chunk by chunk; with ``backward`` the
        gradient is accumulated chunk by chunk (memory bounded by one chunk)."""
        total = 0.0
        for sub in data.chunks(idx):
            c = self._nll(V, mu, data, data.batch(sub, cols=cols), cols) * scale
            if backward:
                c.backward()
            total += float(c.detach())
        return total

    # ---------------------------------------------------------------------------------
    # Initialisation
    # ---------------------------------------------------------------------------------
    def _initial_point(self, data: _Data) -> tuple[torch.Tensor, torch.Tensor]:
        d, p = self.n_pc, data.p
        # saturated parameter of the column means: intercept-only MLE when T(x)=x
        mu0 = self.exponential_family.invert_g(self._colmean.reshape(1, -1)).reshape(-1)
        if isinstance(self.init, (np.ndarray, torch.Tensor)):
            V0 = torch.linalg.qr(torch.as_tensor(self.init).to(self.device, self.dtype))[0]
        elif self.init == "spectral":  # randomized SVD on a memory-capped cell subset
            k = int(min(data.n, 2**13, max(4 * d, 2e8 // p)))
            b = data.batch(np.sort(self._rng.choice(data.n, k, replace=False)), need_X=False)
            Th = data.theta0 + b.delta.to_dense() if b.delta is not None else b.theta
            if b.o is not None:
                Th = Th - b.o[:, None] * (Th != data.theta0)
            V0 = torch.svd_lowrank(Th - mu0, q=min(d + 10, *Th.shape), niter=4)[2][:, :d]
        elif self.init == "random":
            V0 = torch.linalg.qr(torch.randn(p, d, device=self.device, dtype=self.dtype))[0]
        else:
            raise ValueError(f"Unknown init {self.init}")
        return V0.contiguous(), mu0

    # ---------------------------------------------------------------------------------
    # Optimisation
    # ---------------------------------------------------------------------------------
    def _saturated_loading_iter(self, data: _Data) -> tuple[torch.Tensor, torch.Tensor]:
        if self.optimizer == "lbfgs":
            return self._fit_lbfgs(data)
        self.learning_rate_ = self.initial_learning_rate_
        while self.learning_rate_ >= LEARNING_RATE_LIMIT:
            out = self._fit_geoopt(data)
            if out is not None:
                return out
            self._log("\tRESTART BECAUSE INF/NAN FOUND")
            self.learning_rate_ *= self.gamma
        raise ValueError("LEARNING RATE IS TOO SMALL : DID NOT CONVERGE")

    def _fit_geoopt(self, data: _Data):
        import geoopt

        V0, mu0 = self._initial_point(data)
        V = geoopt.ManifoldParameter(V0.clone(), manifold=geoopt.EuclideanStiefel())
        mu = torch.nn.Parameter(mu0.clone())
        lr_V = self.learning_rate_ / np.sqrt(data.p)
        opt = geoopt.optim.RiemannianAdam(
            [{"params": [V], "lr": lr_V},
             {"params": [mu], "lr": self.learning_rate_ * self.intercept_lr_ratio}],
            stabilize=50,
        )
        exact = self.fit_mode == "exact"
        sched = torch.optim.lr_scheduler.ReduceLROnPlateau(
            opt, factor=self.gamma, threshold=self.tol, threshold_mode="abs",
            patience=self.patience * (10 if exact else 1), eps=0.0,
        )
        lr_stop = lr_V * self.gamma ** (self.n_lr_decay + 0.5)
        batch_size = data.n if exact else self.batch_size
        scores, rates = [], []
        self._loadings_epochs, self._intercept_epochs = [], []
        pbar = tqdm(range(self.max_iter), disable=not self.verbose)
        for epoch in pbar:
            epoch_cost = 0.0
            for idx in data.minibatches(batch_size, self._rng):
                opt.zero_grad(set_to_none=True)
                cost = self._accumulate(V, mu, data, idx, self._feature_sample(data.p),
                                        backward=True, scale=1 / len(idx))
                if not np.isfinite(cost) or not torch.isfinite(V).all():
                    return None
                try:  # the Stiefel retraction runs an SVD: a bad step raises here
                    opt.step()
                except torch.linalg.LinAlgError:
                    return None  # -> restart with learning_rate * gamma
                scores.append(cost)
                rates.append(opt.param_groups[0]["lr"])
                epoch_cost += cost * len(idx) / data.n
            pbar.set_postfix(cost=f"{epoch_cost:.4f}")
            if self.store_epochs:
                self._loadings_epochs.append(V.detach().cpu().clone())
                self._intercept_epochs.append(mu.detach().cpu().clone())
            if epoch + 1 >= self.min_epochs:
                sched.step(epoch_cost / data.p)  # mean NLL per matrix entry
                if opt.param_groups[0]["lr"] < lr_stop:
                    break
        self.n_epochs_ = epoch + 1
        self.loadings_learning_scores_.append(scores)
        self.loadings_learning_rates_.append(rates)
        return V.detach(), mu.detach()

    def _fit_lbfgs(self, data: _Data):
        """Full-batch L-BFGS (strong-Wolfe) on the Stiefel manifold via V = qf(W),
        the QR retraction of an unconstrained W (O(p d^2) per call); gradients are
        exact and accumulated chunk by chunk, so memory is bounded by one chunk."""
        V0, mu0 = self._initial_point(data)
        W0 = torch.nn.Parameter(V0.clone())
        mu = torch.nn.Parameter(mu0.clone())

        def stiefel():
            Q, R = torch.linalg.qr(W0)
            return Q * torch.sign(torch.diagonal(R))          # unique: diag(R) > 0

        cols, rows = self._feature_sample(data.p), np.arange(data.n)  # fixed subset
        opt = torch.optim.LBFGS([W0, mu], lr=1, max_iter=self.max_iter, history_size=20,
                                tolerance_grad=1e-7, tolerance_change=1e-12,
                                line_search_fn="strong_wolfe")
        history: list[float] = []

        def closure():
            opt.zero_grad(set_to_none=True)
            W = stiefel()                                      # graph through the map
            V = W.detach().requires_grad_(True)                # chunk grads land here
            cost = self._accumulate(V, mu, data, rows, cols, backward=True, scale=1 / data.n)
            W.backward(V.grad)                                 # one pass through the map
            history.append(cost)
            if np.isfinite(cost) and cost <= best["cost"]:     # keep best finite point
                best.update(cost=cost, V=V.detach().clone(), mu=mu.detach().clone())
            return cost

        best = {"cost": np.inf, "V": V0.clone(), "mu": mu0.clone()}
        try:
            opt.step(closure)
        except RuntimeError as e:  # strong-Wolfe trial step overflowed
            self._log(f"\tL-BFGS stopped early ({e}); returning best point")
        self.loadings_learning_scores_.append(history)
        self.loadings_learning_rates_.append([])
        return best["V"], best["mu"]

    @torch.no_grad()
    def _whiten(self, V, mu, data: _Data):
        r"""Definition 2.8 from the d x d Gram matrix of the scores: the right singular
        vectors of (Theta - 1mu^T) V V^T are V R, R the eigenvectors of U^T U.
        With center_scores, mu <- mu + V m (m = mean score; likelihood unchanged, since mu
        is identified only outside span(V)) and the centred Gram matrix is used.
        Returns (V R, singular values, mu)."""
        G = torch.zeros(self.n_pc, self.n_pc, dtype=torch.float64, device=self.device)
        m = torch.zeros(self.n_pc, dtype=torch.float64, device=self.device)
        for sub in data.chunks(np.arange(data.n)):
            U = data.scores(V, mu, data.batch(sub, need_X=False)).double()
            G += U.T @ U
            m += U.sum(0)
        if self.center_scores:
            m /= data.n
            G -= data.n * torch.outer(m, m)
            mu = mu + V @ m.to(self.dtype)
        evals, R = torch.linalg.eigh(G)
        evals, R = evals.flip(0), R.flip(1)  # descending order
        return V @ R.to(self.dtype), evals.clamp_min(0).sqrt().to(self.dtype), mu
