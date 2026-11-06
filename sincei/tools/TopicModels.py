"""Topic models of a cell-by-region count matrix: LSA and LDA.

`LSA` is the truncated SVD of the TF-IDF equivalent to gensim's
`TfidfModel(smartirs="lfu")`, solved with ARPACK (`scipy.sparse.linalg.svds`).
`LDA` is latent Dirichlet allocation fitted by collapsed Gibbs sampling with
tomotopy.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import anndata as ad
import numpy as np
import tomotopy as tp
from scipy.sparse import csr_matrix, issparse
from scipy.sparse.linalg import svds

from sincei.cli._common_args import available_processors

if TYPE_CHECKING:
    from scipy.sparse import spmatrix

    CountMatrix = ad.AnnData | np.ndarray | spmatrix

DEFAULT_RANDOM_STATE = 42
# Slope of gensim's pivoted unique normalization, `smartirs="lfu"`.
GENSIM_SLOPE = 0.25
LOG_INTERVAL = 10


def count_matrix(X: CountMatrix, binarize: bool = False) -> csr_matrix:
    """The counts of `X`, cells in rows, as CSR. With `binarize`, a non-zero is 1."""
    counts = X.X if isinstance(X, ad.AnnData) else X
    matrix = csr_matrix(counts if issparse(counts) else np.asarray(counts))
    if binarize:
        return csr_matrix(
            ((matrix.data != 0).astype(np.float32), matrix.indices, matrix.indptr),
            shape=matrix.shape,
        )
    return matrix


def idf(holders: np.ndarray, n_cells: int) -> tuple[np.ndarray, float]:
    """The idf of the features, and the pivot of the normalization.

    `holders` are the cells with non-zero counts per feature. The idf is
    `log2(n_cells / holders)`, zero for a feature that no cell holds, and the
    pivot is the mean number of features that a cell holds.
    """
    region_idf = np.where(
        holders > 0, np.log2(n_cells / np.clip(holders, 1.0, None)), 0.0
    )
    return region_idf, float(holders.sum()) / n_cells


def tfidf(matrix: csr_matrix, region_idf: np.ndarray, pivot: float) -> csr_matrix:
    r"""The TF-IDF of the CSR `matrix`, computed on its stored entries.

    It is gensim's `TfidfModel(normalize=True, smartirs="lfu")`, for a count above
    zero:

        (1 + log2(count)) · idf / ((1 - GENSIM_SLOPE) · pivot + GENSIM_SLOPE · kept)

    with `kept` the regions of the cell with a count above zero and an idf above zero.
    A zero stays zero, so the result has the same entries as `matrix`.
    """
    counts = np.asarray(matrix.data, dtype=np.float64)
    rows = np.repeat(np.arange(matrix.shape[0]), np.diff(matrix.indptr))
    held = counts > 0
    weight = np.where(
        held,
        (1.0 + np.log2(np.where(held, counts, 1.0))) * region_idf[matrix.indices],
        0.0,
    )
    kept = np.bincount(
        rows[held & (region_idf[matrix.indices] > 0)], minlength=matrix.shape[0]
    )
    norm = (1.0 - GENSIM_SLOPE) * pivot + GENSIM_SLOPE * kept
    values = weight / norm[rows]
    return csr_matrix(
        (values.astype(np.float32), matrix.indices, matrix.indptr), shape=matrix.shape
    )


def largest_entry_signs(basis: np.ndarray) -> np.ndarray:
    """The sign that makes the entry of largest magnitude of every column positive.

    A zero column keeps its sign. Multiplying the columns by it fixes the arbitrary sign
    of a singular vector.
    """
    largest = np.abs(basis).argmax(axis=0)
    signs = np.sign(basis[largest, np.arange(basis.shape[1])])
    signs[signs == 0] = 1
    return signs


class LSA:
    r"""Latent semantic analysis of a cell-by-region count matrix.

    The counts are turned into the TF-IDF of gensim's `TfidfModel(smartirs="lfu")`
    (`tfidf`), and the components are the top `n_pc` right singular vectors of the
    TF-IDF matrix `W`. The scores of a cell are `w V`.

    Parameters
    ----------
    n_pc : int
        Number of components. Defaults to 20.
    binarize : bool
        Whether to read every value that is not zero as 1. Defaults to False.
    random_state : int or None
        Seed of the start of the ARPACK solver. Defaults to DEFAULT_RANDOM_STATE (42).
    """

    def __init__(
        self,
        n_pc: int = 20,
        binarize: bool = False,
        random_state: int | None = DEFAULT_RANDOM_STATE,
    ) -> None:
        if n_pc < 1:
            msg = f"n_pc={n_pc} is not valid. Use n_pc of 1 or more."
            raise ValueError(msg)
        self.n_pc = n_pc
        self.binarize = binarize
        self.random_state = random_state

        # Set by fit. `V`, regions by components, ordered by singular value, with
        # the entry of largest absolute value of every component positive.
        self._components: np.ndarray | None = None
        # The singular value of every component.
        self._singular_values: np.ndarray | None = None
        # The idf of every region, which transform uses for new cells.
        self._idf: np.ndarray | None = None
        # The pivot of the normalization.
        self._pivot: float | None = None

    def fit(self, X: CountMatrix) -> LSA:
        """Fits the components to `X`, cells in rows and regions in columns."""
        self.fit_transform(X)
        return self

    def fit_transform(self, X: CountMatrix) -> np.ndarray:
        """Fits the components to `X` and returns the scores of its cells."""
        counts = count_matrix(X, self.binarize)
        n, p = counts.shape
        if self.n_pc >= min(n, p):
            msg = (
                f"n_pc={self.n_pc} components need more cells and regions than that, "
                f"but the input has {n} cells and {p} regions. Use n_pc="
                f"{min(n, p) - 1} or fewer."
            )
            raise ValueError(msg)
        if bool((counts.data < 0).any()):
            msg = "LSA takes counts, but the input has negative values."
            raise ValueError(msg)

        holders = np.bincount(counts.indices[counts.data > 0], minlength=p)
        self._idf, self._pivot = idf(holders.astype(np.float64), n)
        weighted = tfidf(counts, self._idf, self._pivot)
        u, s, vt = svds(weighted, k=self.n_pc, rng=self.random_state)
        order = np.argsort(s)[::-1]
        components = vt[order].T.astype(np.float32)
        scores = (u[:, order] * s[order]).astype(np.float32)

        signs = largest_entry_signs(components)
        self._components = components * signs
        self._singular_values = s[order].astype(np.float32)
        return scores * signs

    def transform(self, X: CountMatrix) -> np.ndarray:
        """The scores of the cells of `X`, with the idf of the fit."""
        if self._components is None:
            msg = "LSA is not fitted. Call fit() before transform()."
            raise RuntimeError(msg)
        assert self._idf is not None
        assert self._pivot is not None
        weighted = tfidf(count_matrix(X, self.binarize), self._idf, self._pivot)
        return np.asarray(weighted @ self._components, dtype=np.float32)


class LDA:
    r"""Latent Dirichlet allocation of a cell-by-region count matrix, fitted by
    collapsed Gibbs sampling with tomotopy.

    Each cell is a document and each count of a region is one token of that region.
    Every token holds a topic assignment, and each iteration samples every assignment
    again from the counts of the others (Griffiths and Steyvers, 2004).

    The results are from the last sample: the topic proportions of a cell are
    :math:`(n_{dk} + \alpha) / (n_d + K \alpha)` and the region distribution of a topic
    is :math:`(n_{kw} + \eta) / (n_k + V \eta)`, over every region of the input. A cell
    without counts gets the prior mean, 1/K for every topic.

    Parameters
    ----------
    n_topics : int
        Number of topics K. Defaults to 20.
    alpha : float or None
        Prior of the topic proportions of a cell, the same for every topic. None takes
        1/K. Defaults to None.
    eta : float
        Prior of the region weights of a topic. Defaults to 0.1.
    iterations : int
        Number of Gibbs sampling iterations over all the tokens. Defaults to 1000.
    binarize : bool
        Whether to count every non-zero entry as 1. Defaults to False.
    n_jobs : int or None
        Number of threads for tomotopy. The result of a seed depends on it. None uses
        all available processors. Defaults to None.
    random_state : int
        Seed of the sampler. Defaults to DEFAULT_RANDOM_STATE (42).
    """

    def __init__(
        self,
        n_topics: int = 20,
        alpha: float | None = None,
        eta: float = 0.1,
        iterations: int = 1000,
        binarize: bool = False,
        n_jobs: int | None = None,
        random_state: int = DEFAULT_RANDOM_STATE,
    ) -> None:
        if n_jobs is None:
            n_jobs = available_processors()
        for name, value in (
            ("n_topics", n_topics),
            ("iterations", iterations),
            ("n_jobs", n_jobs),
        ):
            if value < 1:
                msg = f"{name} must be at least 1, not {value}."
                raise ValueError(msg)
        if (alpha is not None and alpha <= 0) or eta <= 0:
            msg = f"alpha={alpha} and eta={eta} are not valid. Use positive values."
            raise ValueError(msg)
        self.n_topics = n_topics
        self.alpha = alpha
        self.eta = eta
        self.iterations = iterations
        self.binarize = binarize
        self.n_jobs = n_jobs
        self.random_state = random_state

        # Set by fit. The region distribution of every topic, topics by regions.
        self._topics: np.ndarray | None = None
        # The log-likelihood per token, every LOG_INTERVAL iterations.
        self._log_likelihood: list[float] = []

    def fit(self, X: CountMatrix) -> LDA:
        """Fits the model to `X`, cells in rows and regions in columns."""
        self.fit_transform(X)
        return self

    def fit_transform(self, X: CountMatrix) -> np.ndarray:
        """Fits the model, and returns the topic proportions of the cells."""
        counts = count_matrix(X, self.binarize)
        if counts.nnz and not np.array_equal(counts.data, np.round(counts.data)):
            msg = (
                "Collapsed Gibbs sampling needs integer counts. Use binarize or round "
                "the matrix."
            )
            raise ValueError(msg)
        if bool((counts.data < 0).any()):
            msg = "LDA takes counts, but the input has negative values."
            raise ValueError(msg)
        alpha = 1.0 / self.n_topics if self.alpha is None else float(self.alpha)
        model = tp.LDAModel(
            tw=tp.TermWeight.ONE,
            min_cf=0,
            min_df=0,
            rm_top=0,
            k=self.n_topics,
            alpha=alpha,
            eta=self.eta,
            seed=self.random_state,
        )
        n_cells, n_regions = counts.shape
        names = np.array([str(region) for region in range(n_regions)])
        values = counts.data.astype(np.int64)
        modelled = []
        for row in range(n_cells):
            low, high = counts.indptr[row], counts.indptr[row + 1]
            tokens = np.repeat(names[counts.indices[low:high]], values[low:high])
            if len(tokens):
                model.add_doc(tokens.tolist())
                modelled.append(row)

        self._log_likelihood = []
        done = 0
        while done < self.iterations:
            step = min(LOG_INTERVAL, self.iterations - done)
            model.train(step, workers=self.n_jobs)
            done += step
            self._log_likelihood.append(float(model.ll_per_word))
            logging.info(
                "Iteration %d: log-likelihood per token %.6f", done, model.ll_per_word
            )

        theta = np.full((n_cells, self.n_topics), 1.0 / self.n_topics, np.float32)
        theta[np.asarray(modelled, dtype=np.int64)] = np.array(
            [document.get_topic_dist() for document in model.docs], dtype=np.float32
        )
        self._topics = self._full_topics(model, n_regions)
        return theta

    def _full_topics(self, model: tp.LDAModel, n_regions: int) -> np.ndarray:
        """(n_kw + eta) / (n_k + V eta) over every region, with or without tokens."""
        regions = np.array([int(word) for word in model.used_vocabs], dtype=np.int64)
        counts = np.zeros((self.n_topics, n_regions))
        for topic in range(self.n_topics):
            raw = np.asarray(model.get_topic_word_dist(topic, normalize=False))
            counts[topic, regions] = raw - self.eta
        weights = np.clip(counts, 0.0, None) + self.eta
        return (weights / weights.sum(axis=1, keepdims=True)).astype(np.float32)
