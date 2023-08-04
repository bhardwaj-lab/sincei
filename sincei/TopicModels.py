# topic models
import numpy as np
import pandas as pd
from gensim import corpora, matutils, models
import copy

# Louvain clustering and UMAP
from networkx import convert_matrix
from sklearn.metrics import pairwise_distances
import leidenalg as la

# import community
import umap
from scanpy._utils import get_igraph_from_adjacency
from scanpy.neighbors import (
    _compute_connectivities_umap,
    _get_indices_distances_from_dense_matrix,
)

### ------ Functions ------


class TOPICMODEL:
    r"""
    Computes LSA or LDA for a given matrix and returns the cell-topic matrix

    Parameters
    ----------
    mat : numpy array
        Matrix of shape (cells, regions)

    cells : list
        List of cells

    regions : list
        List of regions

    nTopics : int
        Number of topics

    smartCode : str
        Smart code for tfidf

    Returns
    -------
    corpus_lsi : gensim corpus
        LSA corpus

    cell_topic : pandas dataframe
        Cell-topic matrix

    corpus_tfidf : gensim corpus
        TFIDF corpus
    """

    def __init__(
        self,
        mat,
        cells,
        regions,
        n_topics,
        smart_code="lfu",
        n_passes=1,
        n_workers=1,
    ):
        self.n_topics = n_topics
        self.smart_code = smart_code
        self.cells = cells
        self.regions_dict = corpora.dictionary.Dictionary([regions])
        self.corpus = matutils.Sparse2Corpus(mat)
        self.n_passes = n_passes
        self.n_workers = n_workers
        self.lsi_model = None
        self.lda_model = None
        self.cell_topic_dist = None
        self.topic_region_dist = None

    def runLSA(self):
        r"Computes LSA for a given matrix and returns the updated object"

        # LSA
        tfidf = models.TfidfModel(self.corpus, id2word=self.regions_dict, normalize=True, smartirs=self.smart_code)
        corpus_tfidf = tfidf[self.corpus]
        self.lsi_model = models.LsiModel(corpus_tfidf, id2word=self.regions_dict, num_topics=self.n_topics)
        self.cell_topic_dist = self.lsi_model[corpus_tfidf]

        # Compute Coherence Score
        coherence_model_lsa = models.CoherenceModel(
            model=self.lsi_model, corpus=self.corpus, dictionary=self.regions_dict, coherence="u_mass"
        )
        coherence_lsa = coherence_model_lsa.get_coherence()
        print("\nCoherence Score: ", coherence_lsa)

    def runLDA(self):
        r"Computes LDA model for a given matrix and returns the updated object"

        self.lda_model = models.LdaMulticore(
            corpus=self.corpus, num_topics=self.n_topics, passes=self.n_passes, workers=self.n_workers
        )
        # get topic distributions for each document
        self.cell_topic_dist = self.lda_model[self.corpus]
        # get topic-word distributions
        self.topic_region_dist = self.lda_model.get_topics()

    def get_cell_topic(self, pop_sparse_cells=False):
        r"Returns cell-topic matrix from the updated object"

        cells = copy.deepcopy(self.cells)
        ## make cell-topic df
        li = [[tup[0] for tup in x] for x in self.cell_topic_dist]
        li_val = [[tup[1] for tup in x] for x in self.cell_topic_dist]

        # if all documents don't have same set of topics, (optionally) remove them
        if len(set([len(x) for x in li_val])) > 1:
            bad_idx = sorted([i for i, v in enumerate(li_val) if len(v) != self.n_topics], reverse=True)
            print("{} Cells were detected which don't contribute to all {} topics.".format(len(bad_idx), self.n_topics))
            if pop_sparse_cells:
                print("Removing these cells from the analysis")
                [li_val.pop(x) for x in bad_idx]
                [li.pop(x) for x in bad_idx]
                [cells.pop(x) for x in bad_idx]
            else:
                print("Not implemented! Need to fill these entries with zeros")

        li_val = np.stack(li_val)
        cell_topic = pd.DataFrame(li_val, columns=li[0])
        cell_topic.index = cells

        return cell_topic
