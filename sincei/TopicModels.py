# topic models
import numpy as np
import pandas as pd
from gensim import corpora, matutils, models
import copy

### ------ Functions ------

class TOPICMODEL:
    r"""
    Computes LSA or LDA for a given matrix and returns the cell-topic matrix

    Parameters
    ----------
    mat : scipy.sparse.csr_matrix
        Sparse matrix of shape (cells, regions)

    cells : list
        List of Cell IDs (corresponding to the input matrix rows)

    regions : list
        List of Regions (corresponding to the input matrix columns)

    n_topics : int
        Number of Topics / Principal Components

    smart_code : str
        SMART (System for the Mechanical Analysis and Retrieval of Text) code for weighting of input matrix for TFIDF.
        Only valid for the LSA model. The default ("lfu") corresponds to "log"TF * IDF, and "pivoted unique" normalization of document length. For more information, see: https://en.wikipedia.org/wiki/SMART_Information_Retrieval_System

    Returns
    -------
    An object of class TOPICMODELS.
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
        self.shape = (len(cells), len(regions))

    def runLSA(self):
        r"""
        Computes LSA for a given matrix and returns the updated object

        Returns
        -------
        corpus_tfidf : gensim corpus
            TFIDF normalized corpus
        corpus_lsi : gensim corpus
            LSA corpus
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
        Computes LDA model for a given matrix and returns the updated object

        Returns
        -------
        cell_topic : pandas dataframe
            Cell-topic matrix
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
        Get cell-topic matrix from the updated object

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
            print("{} Cells were detected which don't contribute to all {} topics.".format(len(bad_idx), self.n_topics))
            if pop_sparse_cells:
                print("Removing these cells from the analysis")
                for x in bad_idx:
                    li_val.pop(x)
                    li.pop(x)
                    cells.pop(x)
            else:
                print("Not implemented! Need to fill these entries with zeros")

        li_val = np.stack(li_val)
        cell_topic = pd.DataFrame(li_val, columns=li[0])
        cell_topic.index = cells

        return cell_topic
