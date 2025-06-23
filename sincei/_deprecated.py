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
