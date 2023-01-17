# topic models
import numpy as np
import pandas as pd
import scanpy as sc
from gensim import corpora, matutils, models

# Louvain clustering and UMAP
from networkx import convert_matrix
from sklearn.metrics import pairwise_distances
import leidenalg as la
import community
import umap
from scanpy._utils import get_igraph_from_adjacency
from scanpy.neighbors import _compute_connectivities_umap,  _get_indices_distances_from_dense_matrix

### ------ Functions ------


def LSA_gensim(mat, cells, regions, nTopics, smartCode='lfu'):
    '''
    Perform latent schemantic analysis on a gene*cell sparse matrix (genes=rows, cells=columns)
    Return:
        corpus_lsi: trained LSI model
        cell_topic: cell*topic matrix
        corpus_tfidf: TF-IDF transformed matrix
    '''
    # LSA
    regions_dict = corpora.dictionary.Dictionary([regions])
    corpus = matutils.Sparse2Corpus(mat)
    tfidf = models.TfidfModel(corpus, id2word=regions_dict, normalize=True, smartirs=smartCode)
    corpus_tfidf = tfidf[corpus]
    lsi_model = models.LsiModel(corpus_tfidf, id2word=regions_dict, num_topics=nTopics)
    corpus_lsi = lsi_model[corpus_tfidf]

    # Compute Coherence Score
    coherence_model_lsa = models.CoherenceModel(model=lsi_model, corpus=corpus,
                                                   dictionary=regions_dict, coherence='u_mass')
    coherence_lsa = coherence_model_lsa.get_coherence()
    print('\nCoherence Score: ', coherence_lsa)

    ## make cell-topic df
    li = [[tup[0] for tup in x] for x in corpus_lsi]
    li_val = [[tup[1] for tup in x] for x in corpus_lsi]
    if len(set([len(x) for x in li_val])) > 1: # if all documents don't have same set of topics
        bad_idx = sorted([i for i,v in enumerate(li_val) if len(v) != nTopics], reverse=True)
        print("{} Cells were detected which don't contribute to all {} topics. Removing them!".format(len(bad_idx), nTopics))
        [li_val.pop(x) for x in bad_idx]
        [li.pop(x) for x in bad_idx]
        [cells.pop(x) for x in bad_idx]
    li_val = np.stack(li_val)
    cell_topic = pd.DataFrame(li_val, columns=li[0])
    cell_topic.index = cells

    return corpus_lsi, cell_topic, corpus_tfidf
