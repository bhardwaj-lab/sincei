import umap
import pandas as pd
import scanpy as sc
from scanpy.neighbors import (
    _compute_connectivities_umap,
    _get_indices_distances_from_dense_matrix,
)
from scanpy._utils import get_igraph_from_adjacency
import igraph as ig
import leidenalg as la
from scipy import sparse

import sys

# own modules
from sincei.TopicModels import TOPICMODEL
from sincei._deprecated import cluster_LSA

# umap.__version__ : should be >= 0.5.1

## literature on multi-graph clustering
"""
https://github.com/crisbodnar/regularised-spectral-clustering
https://arxiv.org/pdf/2103.16534.pdf
https://arxiv.org/pdf/2010.12301.pdf
https://arxiv.org/pdf/2010.15456.pdf
https://chenannie45.github.io/SDM18_MANE.pdf
https://www.researchgate.net/post/Python_or_R_packages_for_multilayer_network
http://www.mkivela.com/pymnet/visualizing.html#visualization-tutorial
https://cran.r-project.org/web/packages/multinet/multinet.pdf
https://sciendo.com/pdf/10.2478/amcs-2019-0010
https://web.media.mit.edu/~xdong/paper/tsp14.pdf
"""


def multiModal_clustering(mode1_adata, mode2_adata, column_key="barcode_nla", nK=20):
    r"""
    Performs multi-graph clustering on matched keys(barcodes) bw two anndata objects.

    Parameters
    ----------
    mode1_adata : AnnData
        AnnData object for mode 1
    mode2_adata : AnnData
        AnnData object for mode 2
    column_key : str
        Column name for the barcode in mode 1 and 2
    nK : int
        Number of clusters to use for clustering

    Returns
    -------
    multi_umap : DataFrame
        DataFrame with UMAP coordinates and cluster labels
    mode1_adata : AnnData
        AnnData object for mode 1
    mode2_adata : AnnData
        AnnData object for mode 2
    """

    # subset RNA anndata
    mode1_adata = mode1_adata[[i for i, v in enumerate(mode1_adata.obs[column_key]) if v in mode2_adata.obs.index]]
    # re-compute distances/clusters/umap
    sc.pp.neighbors(mode1_adata, n_neighbors=nK)
    sc.tl.louvain(mode1_adata)
    sc.tl.umap(mode1_adata, min_dist=0.5, spread=5, init_pos="random", random_state=42)
    # get graph
    G_rna = get_igraph_from_adjacency(mode1_adata.obsp["connectivities"], directed=True)

    # subset chic anndata
    # mode2_adata=mode2_adata[[i for i,v in enumerate(mode2_adata.obs_names.to_list()) if v in set(mode1_adata.obs[column_key])]]
    mode2_adata = mode2_adata[mode1_adata.obs[column_key].tolist()]
    # re-do LSA
    mtx = sparse.csr_matrix(mode2_adata.X.transpose())
    dat = TOPICMODEL(
        mtx,
        list(mode2_adata.obs.index),
        list(mode2_adata.var.index),
        nTopics=20,
        smartCode="lfu",
    )
    dat.runLSA()
    cell_topic = dat.get_cell_topic()

    # get graph
    chic_umap, G_chic = cluster_LSA(cell_topic, modularityAlg="leiden", resolution=1, nk=nK)
    mode2_adata.obsm["X_pca"] = cell_topic
    mode2_adata.obsm["X_umap"] = chic_umap
    # leiden multi-layer clustering
    optimiser = la.Optimiser()
    part_rna = la.ModularityVertexPartition(G_rna)
    part_chic = la.ModularityVertexPartition(G_chic)
    #
    # part_rna = la.CPMVertexPartition(G_rna, resolution=res_layer1)
    # part_chic = la.CPMVertexPartition(G_chic, resolution=res_layer2)
    optimiser.optimise_partition_multiplex([part_rna, part_chic], layer_weights=[2, 1], n_iterations=-1)
    print("Detected clusters: ", set(part_chic.membership))
    chic_umap["cluster_multi"] = part_chic.membership

    # merge RNA and ChIC UMAPs
    rna_umap = pd.DataFrame(mode1_adata.obsm["X_umap"], columns=["RNA_UMAP1", "RNA_UMAP2"])
    rna_umap.index = mode1_adata.obs.index
    rna_umap["cluster_mode1"] = mode1_adata.obs.louvain
    rna_umap[column_key] = mode1_adata.obs[column_key]
    multi_umap = rna_umap.merge(chic_umap, left_index=False, right_index=True, left_on=column_key)
    multi_umap["cluster_mode1"] = [int(x) for x in multi_umap["cluster_mode1"].to_list()]

    return multi_umap, mode1_adata, mode2_adata


## run aligned UMAPs between two dfs with PCs
def umap_aligned(pca_mode1, pca_mode2, nK=15, distance_metric="eucledian"):
    r"""
    Aligns two UMAP embeddings using the UMAP AlignedUMAP class

    Parameters
    ----------
    pca_mode1 : pandas.DataFrame
        UMAP embedding of RNA data
    pca_mode2 : pandas.DataFrame
        UMAP embedding of CHiC data
    nK : int
        Number of nearest neighbors to use for UMAP
    distance_metric : str
        Distance metric to use for UMAP

    Returns
    -------
    pandas.DataFrame
        Aligned UMAP embedding of RNA data
    pandas.DataFrame
        Aligned UMAP embedding of CHiC data

    Examples
    --------

    >>> test = Tester()
    >>> pca_mode1 = test.pca_mode1
    >>> pca_mode2 = test.pca_mode2
    >>> umap_aligned(pca_mode1, pca_mode2)
    (                 aligned_UMAP1  aligned_UMAP2
    0  -0.012995  0.001206
    1   0.
    """

    pca_mode1 = pca_mode1.loc[pca_mode2.index]
    keys = [x for x in range(len(pca_mode1.index))]
    values = []
    for x in pca_mode1.index:
        values.append([i for i, v in enumerate(pca_mode2.index) if v == x])
    values = [x[0] for x in values]
    relation_dict = {i: v for i, v in zip(keys, values)}

    chic_idx = pca_mode2.index
    rna_idx = pca_mode1.index
    pca_mode2 = pca_mode2.reset_index().drop("index", axis=1)
    pca_mode1 = pca_mode1.reset_index().drop("index", axis=1)
    # UMAP
    UA = umap.AlignedUMAP(
        spread=10,
        min_dist=0.01,
        n_neighbors=nK,
        metric=distance_metric,
        init="random",
        random_state=42,
    )
    aligned_umap = UA.fit([pca_mode1, pca_mode2], relations=[relation_dict])

    mode1_umap_aligned = pd.DataFrame(aligned_umap.embeddings_[0])
    mode2_umap_aligned = pd.DataFrame(aligned_umap.embeddings_[1])
    mode2_umap_aligned.index = chic_idx
    mode1_umap_aligned.index = rna_idx
    mode2_umap_aligned.columns = ["aligned_UMAP1", "aligned_UMAP2"]
    mode1_umap_aligned.columns = ["aligned_UMAP1", "aligned_UMAP2"]

    return mode1_umap_aligned, mode2_umap_aligned
