def read_mtx(prefix):
    r"""Reads a matrix from a mtx file

    Parameters
    ----------
    prefix : str
        prefix of the mtx file

    Returns
    -------
    mtx : scipy.sparse.csr_matrix
        matrix
    rownames : list
        list of row names
    colnames : list
        list of column names

    Examples
    --------

    >>> mtx, rownames, colnames = read_mtx("test.mtx")
    """
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
    r"""Preprocesses a sparse matrix for use with scanpy

    Parameters
    ----------
    sparse_mtx : scipy.sparse.csr_matrix
        Sparse matrix of read counts

    rownames : list
        List of row names

    colnames : list
        List of column names

    min_cell_sum : int
        Minimum number of reads per cell

    min_region_sum : int
        Minimum number of reads per region

    Returns
    -------
    scanpy.AnnData
        AnnData object with sparse matrix of read counts


    Examples
    --------

    >>> test = Tester()
    >>> sparse_mtx = test.sparse_mtx
    >>> rownames = test.rownames
    >>> colnames = test.colnames
    >>> adata = preprocess_mtx(sparse_mtx, rownames, colnames, min_cell_sum=1, min_region_sum=1)
    """

    from itertools import compress

    ## binarize
    nonzero_mask = np.array(sparse_mtx[sparse_mtx.nonzero()] > 1)[0]
    rows = sparse_mtx.nonzero()[0][nonzero_mask]
    cols = sparse_mtx.nonzero()[1][nonzero_mask]
    sparse_mtx[rows, cols] = 1

    ## filter low counts
    colmask = np.array(np.sum(sparse_mtx, axis=0) >= min_cell_sum)[0]
    rowmask = np.array(np.sum(sparse_mtx, axis=1) >= min_region_sum)
    rowmask = np.array([x[0] for x in rowmask])
    sparse_mtx = sparse_mtx[rowmask, :]
    sparse_mtx = sparse_mtx[:, colmask]

    ## create anndata
    row_subset = list(compress(rownames, rowmask))
    col_subset = list(compress(colnames, colmask))
    adata = anndata.AnnData(
        sparse_mtx.transpose(),
        obs=pd.DataFrame(col_subset),
        var=pd.DataFrame(row_subset),
    )

    return adata


# Imports for cluster_topic
# from sklearn.metrics import pairwise_distances
# from scanpy.neighbors import (
#     _compute_connectivities_umap,
#     _get_indices_distances_from_dense_matrix,
# )
# from scanpy._utils import get_igraph_from_adjacency
# import leidenalg as la
# import umap
# import pandas as pd


def cluster_topic(
    cell_topic,
    distance_metric="cosine",
    nk=30,
    resolution=1.0,
    connectivity_graph=True,
):
    r"""Cluster cells using the output of LSA_gensim or LDA_gensim using the Leiden algorithm.

    Parameters
    ----------
    cell_topic : pandas dataframe
        cell_topic dataframe with cell names as index and topic proportions as columns.

    distance_metric : str
        Distance metric to use. Default: 'cosine'

    nk : int
        Number of nearest neighbors to use. Default: 30

    resolution : float
        Resolution parameter for modularity algorithm. Default: 1.0

    connectivity_graph : bool
        Whether to use a connectivity graph or a distance graph. Default: True

    Returns
    -------
    umap_df : pandas dataframe
        UMAP embedding with cluster labels.

    G : igraph object
        Graph object.
    """

    # cluster on cell-topic dist
    _distances = pairwise_distances(cell_topic.iloc[:, 1:], metric=distance_metric)
    knn_indices, knn_distances = _get_indices_distances_from_dense_matrix(_distances, nk)
    distances, connectivities = _compute_connectivities_umap(knn_indices, knn_distances, _distances.shape[0], nk)

    if connectivity_graph:
        G = get_igraph_from_adjacency(connectivities, directed=True)
    else:
        G = get_igraph_from_adjacency(distances, directed=True)
    partition = la.find_partition(
        G,
        la.RBConfigurationVertexPartition,
        weights="weight",
        seed=42,
        resolution_parameter=resolution,
    )
    cell_topic["cluster"] = partition.membership

    # umap on cell-topic dist
    um = umap.UMAP(
        spread=5,
        min_dist=0.1,
        n_neighbors=nk,
        metric=distance_metric,
        init="random",
        random_state=42,
    )
    umfit = um.fit(cell_topic.iloc[:, 0 : (len(cell_topic.columns) - 1)])
    umap_df = pd.DataFrame(umfit.embedding_)
    umap_df.columns = ["UMAP1", "UMAP2"]
    umap_df["cluster"] = list(cell_topic.cluster)
    umap_df.index = cell_topic.index

    return umap_df, G


## LSA using scipy (older version)
# def runLSA(sparse_mtx, nPCs, scaleFactor):
# tf = (sparse_mtx.transpose() / sparse_mtx.sum(axis=0)).transpose()
#    tf = sparse_mtx / sparse_mtx.sum(axis=0)
#    tf = np.log1p(tf * scaleFactor)

#    idf = np.log1p(sparse_mtx.shape[1] / sparse_mtx.sum(axis=1))
#    tfidf = np.multiply(tf, idf)
#    svd = TruncatedSVD(n_components=nPCs)
#    pca = svd.fit(np.nan_to_num(tfidf, nan = 0.0))

#    return pca, tfidf

## for anndata
# def lsa_anndata(adata, n_pcs, scale_factor):
#   from scipy.sparse import issparse, coo_matrix, csr_matrix
#    mtx = sparse.csr_matrix(adata.X)
#    lsa_out, tfidf = runLSA(mtx.transpose(), n_pcs, scale_factor)

#    adata.X = np.squeeze(np.asarray(tfidf.transpose()))
#    adata.obsm['X_pca'] = lsa_out.components_.transpose()
#    adata.uns['pca_variance'] = lsa_out.explained_variance_
#    adata.uns['pca_variance_ratio'] = lsa_out.explained_variance_ratio_
#    adata.layers['raw_counts'] = mtx

#    return adata

# def UMAP_clustering(adata):

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
