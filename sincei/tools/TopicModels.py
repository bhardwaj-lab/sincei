# topic models
from __future__ import annotations

from typing import TYPE_CHECKING, cast

import numpy as np
import pandas as pd
from gensim import corpora, matutils, models
from sklearn.preprocessing import binarize as sklearn_binarize

if TYPE_CHECKING:
    from collections.abc import Iterable

    import anndata as ad
    from scipy import sparse

    # The matrix layouts `anndata.AnnData.X` holds for a count matrix.
    CountMatrix = (
        np.ndarray
        | sparse.csr_matrix
        | sparse.csc_matrix
        | sparse.csr_array
        | sparse.csc_array
    )

### ------ Functions ------


class TOPICMODEL:
    r"""
    Computes LSA or LDA for a given matrix and returns the cell-topic matrix.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing the data matrix in adata.X, with cells in
        adata.obs_names and regions in adata.var_names.
    n_topics : int
        Number of Topics / Principal Components for modeling.
    binarize : bool, optional
        If True, the input matrix will be binarized (default is False). Recommended
        for LDA.
    smart_code : str
        SMART (System for the Mechanical Analysis and Retrieval of Text) code for
        weighting of input matrix for TFIDF. Only valid for the LSA model. The default
        ("lfu") corresponds to "log" TF * IDF, and "pivoted unique" normalization of
        document length. For more information, see:
        https://en.wikipedia.org/wiki/SMART_Information_Retrieval_System
    n_passes : int, optional
        Number of passes for the LDA model. Default is 1.
    n_workers : int, optional
        Number of workers for the LDA model. Default is 1.
    """

    def __init__(
        self,
        adata: ad.AnnData,
        n_topics: int,
        binarize: bool = False,
        smart_code: str = "lfu",
        n_passes: int = 1,
        n_workers: int = 1,
    ) -> None:
        self.cells = adata.obs_names.to_list()
        self.regions_dict = corpora.dictionary.Dictionary([adata.var_names.to_list()])
        mtx = cast("CountMatrix", adata.X).copy().transpose()
        if binarize:
            mtx = sklearn_binarize(mtx, copy=True)
        self.corpus = matutils.Sparse2Corpus(mtx)
        self.shape = adata.shape
        self.n_topics = n_topics
        self.smart_code = smart_code
        self.n_passes = n_passes
        self.n_workers = n_workers
        self.corpus_tfidf = None
        self.lsi_model: models.LsiModel | None = None
        self.lda_model: models.LdaMulticore | None = None
        self.cell_topic_dist: Iterable[list[tuple[int, float]]] | None = None
        self.topic_region_dist: np.ndarray | None = None

    def runLSA(self) -> None:
        r"""
        Computes LSA for a given matrix and updates the ``TOPICMODEL`` object.
        """

        # LSA
        tfidf = models.TfidfModel(
            self.corpus,
            id2word=self.regions_dict,
            normalize=True,
            smartirs=self.smart_code,
        )
        self.corpus_tfidf = tfidf[self.corpus]
        self.lsi_model = models.LsiModel(
            self.corpus_tfidf, id2word=self.regions_dict, num_topics=self.n_topics
        )
        # lsi[X] computes U^-1*X, which equals V*S (its shape is num_docs * num_topics).
        self.cell_topic_dist = self.lsi_model[self.corpus_tfidf]

        # Compute Coherence Score
        coherence_model_lsa = models.CoherenceModel(
            model=self.lsi_model,
            corpus=self.corpus,
            dictionary=self.regions_dict,
            coherence="u_mass",
        )
        coherence_lsa = coherence_model_lsa.get_coherence()
        print("\nCoherence Score: ", coherence_lsa)

    def runLDA(
        self,
        iterations: int = 500,
        alpha: float = 50.0,
        eta: float = 0.1,
        gamma_threshold: float = 0.001,
    ) -> None:
        r"""
        Computes LDA model for a given matrix and updates the ``TOPICMODEL`` object.

        Parameters
        ----------
        iterations : int, optional
            Maximum number of iterations through the corpus when inferring the topic
            distribution of a corpus. Default is 500.
        alpha : float, optional
            A-priori belief on the cell-topic distribution. Default is 50.0.
        eta : float, optional
            A-priori belief on the topic-region distribution. Default is 0.1.
        gamma_threshold : float, optional
            Minimum change in the value of the gamma parameters to continue iterating.
            Default is 0.001.
        """

        self.lda_model = models.LdaMulticore(
            corpus=self.corpus,
            num_topics=self.n_topics,
            passes=self.n_passes,
            iterations=iterations,
            alpha=alpha,
            eta=eta,
            gamma_threshold=gamma_threshold,
            chunksize=5000,
            workers=self.n_workers,
            minimum_probability=0.0,
        )
        # get topic distributions for each document as dense topic vectors
        self.cell_topic_dist = self.lda_model.get_document_topics(
            self.corpus, minimum_probability=0.0
        )
        # get topic-word distributions
        self.topic_region_dist = self.lda_model.get_topics()

    def get_cell_topic(self) -> pd.DataFrame:
        r"""
        Get cell-topic matrix from the ``TOPICMODEL`` object.

        Returns
        -------
        cell_topic : pandas.DataFrame
            Cell-topic matrix (cells x topics).
        """
        topic_dist = self.cell_topic_dist
        if topic_dist is None:
            msg = "No topic model fitted yet. Call runLSA() or runLDA() first."
            raise RuntimeError(msg)

        weights = np.zeros((len(self.cells), self.n_topics), dtype=float)
        for i, topic_weights in enumerate(topic_dist):
            for topic_id, weight in topic_weights:
                if topic_id < self.n_topics:
                    weights[i, topic_id] = weight

        cell_topic = pd.DataFrame(
            weights, columns=pd.Index([f"topic_{x}" for x in range(self.n_topics)])
        )
        cell_topic.index = pd.Index(self.cells)

        return cell_topic
