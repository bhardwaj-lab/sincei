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

def read_mtx(prefix):
    mtx = io.mmread(prefix + ".counts.mtx")
    mtx = mtx.tocsr()
    with open(prefix + ".colnames.txt") as col:
        colnames = col.read().splitlines()
    col.close()
    with open(prefix + ".rownames.txt") as row:
        rownames = row.read().splitlines()
    row.close()

    return mtx, rownames, colnames

# from mtx
def preprocess_mtx(sparse_mtx, rownames, colnames, min_cell_sum, min_region_sum):
    from itertools import compress
    ## binarize
    nonzero_mask = np.array(sparse_mtx[sparse_mtx.nonzero()] > 1)[0]
    rows = sparse_mtx.nonzero()[0][nonzero_mask]
    cols = sparse_mtx.nonzero()[1][nonzero_mask]
    sparse_mtx[rows, cols] = 1

    ## filter low counts
    colmask = np.array(np.sum(sparse_mtx, axis = 0) >= min_cell_sum)[0]
    rowmask = np.array(np.sum(sparse_mtx, axis = 1) >= min_region_sum)
    rowmask = np.array([x[0] for x in rowmask])
    sparse_mtx = sparse_mtx[rowmask, :]
    sparse_mtx = sparse_mtx[:, colmask]

    ## create anndata
    row_subset = list(compress(rownames, rowmask))
    col_subset = list(compress(colnames, colmask))
    adata = anndata.AnnData(sparse_mtx.transpose(),
                            obs=pd.DataFrame(col_subset),
                            var=pd.DataFrame(row_subset))

    return adata

# from anndata
def preprocess_adata(adata, min_cell_sum, min_region_sum):
    # binarize
    sparse_mtx = adata.X
    nonzero_mask = np.array(sparse_mtx[sparse_mtx.nonzero()] > 1)[0]
    rows = sparse_mtx.nonzero()[0][nonzero_mask]
    cols = sparse_mtx.nonzero()[1][nonzero_mask]
    sparse_mtx[rows, cols] = 1
    adata.x = sparse_mtx
    # save some QC
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    # filter
    adata = adata[adata.obs.n_genes_by_counts >= min_cell_sum, :]
    adata = adata[:, adata.var.n_cells_by_counts >= min_region_sum]

    return adata

def LSA_gensim(mat, cells, regions, nTopics, smartCode='lfu'):
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

def cluster_LSA(cell_topic, modularityAlg = 'leiden', distance_metric='cosine', nk=30, resolution=1.0, connectivity_graph=True):

    # cluster on cel-topic dist
    _distances = pairwise_distances(cell_topic.iloc[:, 1:], metric=distance_metric)
    knn_indices, knn_distances = _get_indices_distances_from_dense_matrix(_distances, nk)
    distances, connectivities = _compute_connectivities_umap(knn_indices,
                                                             knn_distances,
                                                             _distances.shape[0], nk)


    if modularityAlg == 'leiden':
        if connectivity_graph:
            G = get_igraph_from_adjacency(connectivities, directed=True)
        else:
            G = get_igraph_from_adjacency(distances, directed=True)
        partition = la.find_partition(G,
                              la.RBConfigurationVertexPartition,
                              weights='weight',
                              seed=42,
                              resolution_parameter=resolution)
        cell_topic['cluster'] = partition.membership
    else:
        if connectivity_graph:
            G = convert_matrix.from_numpy_array(connectivities)
        else:
            G = convert_matrix.from_numpy_array(distances)
        partition = community.best_partition(G, resolution=resolution, random_state=42)
        cell_topic['cluster'] = partition.values()

    # umap on cell-topic dist
    um = umap.UMAP(spread = 5, min_dist=0.1, n_neighbors=nk, metric=distance_metric, init='random', random_state=42)
    umfit = um.fit(cell_topic.iloc[:, 0:(len(cell_topic.columns) - 1)])
    umap_df = pd.DataFrame(umfit.embedding_)
    umap_df.columns = ['UMAP1', 'UMAP2']
    umap_df['cluster'] = list(cell_topic.cluster)
    umap_df.index = cell_topic.index

    return umap_df, G

## LSA using scipy (older version)
#def runLSA(sparse_mtx, nPCs, scaleFactor):
    #tf = (sparse_mtx.transpose() / sparse_mtx.sum(axis=0)).transpose()
#    tf = sparse_mtx / sparse_mtx.sum(axis=0)
#    tf = np.log1p(tf * scaleFactor)

#    idf = np.log1p(sparse_mtx.shape[1] / sparse_mtx.sum(axis=1))
#    tfidf = np.multiply(tf, idf)
#    svd = TruncatedSVD(n_components=nPCs)
#    pca = svd.fit(np.nan_to_num(tfidf, nan = 0.0))

#    return pca, tfidf

## for anndata
#def lsa_anndata(adata, n_pcs, scale_factor):
#   from scipy.sparse import issparse, coo_matrix, csr_matrix
#    mtx = sparse.csr_matrix(adata.X)
#    lsa_out, tfidf = runLSA(mtx.transpose(), n_pcs, scale_factor)

#    adata.X = np.squeeze(np.asarray(tfidf.transpose()))
#    adata.obsm['X_pca'] = lsa_out.components_.transpose()
#    adata.uns['pca_variance'] = lsa_out.explained_variance_
#    adata.uns['pca_variance_ratio'] = lsa_out.explained_variance_ratio_
#    adata.layers['raw_counts'] = mtx

#    return adata

#def UMAP_clustering(adata):

    # make umap on PCA
#    lsa_out = adata.obsm['X_pca'].transpose()[1:, :]
#    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, spread=1.0, metric='euclidean', init = 'random')
#    embeddings = reducer.fit_transform(lsa_out.transpose())
#    adata.obsm['X_umap'] = embeddings

    # louvain (from scanpy)
#    scp.pp.neighbors(adata)
#    scp.tl.louvain(adata)
#    cluster_id = [int(x) for x in adata.obs['louvain'].to_list()]

    # return
#    return adata
