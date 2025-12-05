# topic models
import numpy as np
import pandas as pd
from gensim import corpora, matutils, models
import copy
from sklearn.preprocessing import binarize

### ------ Functions ------


class TOPICMODEL:
    r"""
    Computes LSA or LDA for a given matrix and returns the cell-topic matrix.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing the data matrix in adata.X, with cells in adata.obs_names and regions in adata.var_names.
    n_topics : int
        Number of Topics / Principal Components for modeling.
    binarize : bool, optional
        If True, the input matrix will be binarized (default is False). Recommended for LDA.
    smart_code : str
        SMART (System for the Mechanical Analysis and Retrieval of Text) code for weighting of input matrix for TFIDF.
        Only valid for the LSA model. The default ("lfu") corresponds to "log"TF * IDF, and "pivoted unique" normalization of document length. For more information, see: https://en.wikipedia.org/wiki/SMART_Information_Retrieval_System
    n_passes : int, optional
        Number of passes for the LDA model. Default is 1.
    n_workers : int, optional
        Number of workers for the LDA model. Default is 1.
    """

    def __init__(
        self,
        adata,
        n_topics,
        binarize=False,
        smart_code="lfu",
        n_passes=1,
        n_workers=1,
    ):
        self.cells = adata.obs_names.to_list()
        self.regions_dict = corpora.dictionary.Dictionary([adata.var_names.to_list()])
        mtx = adata.X.copy().transpose()
        if binarize:
            mtx = binarize(mtx, copy=True)
        self.corpus = matutils.Sparse2Corpus(mtx)
        self.shape = adata.shape
        self.n_topics = n_topics
        self.smart_code = smart_code
        self.n_passes = n_passes
        self.n_workers = n_workers
        self.lsi_model = None
        self.lda_model = None
        self.cell_topic_dist = None
        self.topic_region_dist = None

    def runLSA(self):
        r"""
        Computes LSA for a given matrix and updates the ``TOPICMODEL`` object.
        """

        # LSA
        tfidf = models.TfidfModel(self.corpus, id2word=self.regions_dict, normalize=True, smartirs=self.smart_code)
        self.corpus_tfidf = tfidf[self.corpus]
        self.lsi_model = models.LsiModel(self.corpus_tfidf, id2word=self.regions_dict, num_topics=self.n_topics)
        self.cell_topic_dist = self.lsi_model[
            self.corpus_tfidf
        ]  # lsi[X] computes U^-1*X, which equals V*S (its shape is num_docs * num_topics).

        # Compute Coherence Score
        coherence_model_lsa = models.CoherenceModel(
            model=self.lsi_model, corpus=self.corpus, dictionary=self.regions_dict, coherence="u_mass"
        )
        coherence_lsa = coherence_model_lsa.get_coherence()
        print("\nCoherence Score: ", coherence_lsa)

    def runLDA(self):
        r"""
        Computes LDA model for a given matrix and updates the ``TOPICMODEL`` object.
        """

        self.lda_model = models.LdaMulticore(
            corpus=self.corpus, num_topics=self.n_topics, passes=self.n_passes, workers=self.n_workers
        )
        # get topic distributions for each document
        self.cell_topic_dist = self.lda_model[self.corpus]
        # get topic-word distributions
        self.topic_region_dist = self.lda_model.get_topics()

    def get_cell_topic(self, pop_sparse_cells=False):
        r"""
        Get cell-topic matrix from the ``TOPICMODEL`` object.

        Returns
        -------
        cell_topic : pandas dataframe
            Cell-topic matrix
        """

        cells = copy.deepcopy(self.cells)
        ## make cell-topic df
        li = [[tup[0] for tup in x] for x in self.cell_topic_dist]
        li_val = [[tup[1] for tup in x] for x in self.cell_topic_dist]

        # if all documents don't have same set of topics, (optionally) remove them
        if len(set([len(x) for x in li_val])) > 1:
            bad_idx = sorted([i for i, v in enumerate(li_val) if len(v) != self.n_topics], reverse=True)
            print(f"{len(bad_idx)} cells were detected which don't contribute to all {self.n_topics} topics.")
            if pop_sparse_cells:
                print("Removing these cells from the analysis")
                for x in bad_idx:
                    li_val.pop(x)
                    li.pop(x)
                    cells.pop(x)

                li_val = np.stack(li_val)
                cell_topic = pd.DataFrame(li_val, columns=[f"topic_{x}" for x in range(self.n_topics)])
            else:
                cell_topic = np.zeros((len(li_val), self.n_topics))
                for i, v in enumerate(li_val):
                    for j, val in enumerate(v):
                        print(f"Index [{i}, {li[i][j]}] = {val}")
                        cell_topic[i, li[i][j]] = val
                cell_topic = pd.DataFrame(cell_topic, columns=[f"topic_{x}" for x in range(self.n_topics)])

        else:
            li_val = np.stack(li_val)
            cell_topic = pd.DataFrame(li_val, columns=[f"topic_{x}" for x in range(self.n_topics)])

        cell_topic.index = cells

        return cell_topic
