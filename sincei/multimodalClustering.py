import umap
import pandas as pd
import scanpy as sc

from scanpy._utils import get_igraph_from_adjacency
import leidenalg as la

# own modules
from sincei.TopicModels import TOPICMODEL
from sincei.Utilities import cluster_topic

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


def multiModal_clustering(mdata, modalities=None, topic=None, modal_weights=None, column_key=None, nK=20, n_topics=20):
    r"""
    Performs multi-graph clustering on matched keys(barcodes) of a mudata object.

    Parameters
    ----------
    mudata : MuData
        MuData object containing several data modalities
    modalities : list[str]
        List of modalities to use for clustering, e.g. ["RNA", "CHiC"]
    topic : list[str]
        Whether to use topic modeling for each modality. Choose between "False", "LSA" or "LDA". If None, no topic modeling is performed.
    modal_weights : list[float]
        Weights for each modality in the clustering process. Default is equal weighting. E.g. for RNA and CHiC, use [2, 1].
    column_key : str, optional
        Column name for the barcode. If None, the index of obs for each modality is used.
    nK : int
        Number of clusters to use for clustering.
    n_topics : int
        Number of topics to use for topic model. Default is 20.

    Returns
    -------
    multi_umap : DataFrame
        DataFrame with UMAP coordinates and cluster labels
    mudata : MuData
        MuData object containing several data modalities
    """

    # check if modalities are provided, otherwise use all
    if modalities is None:
        raise ValueError(f"Provide modalities to use for clustering.")
    # check if modalities are in mudata object
    for mod in modalities:
        if mod not in mdata.mod.keys():
            raise ValueError(f"Modality {mod} not found in MuData object.")
    # check if topic is provided, otherwise use False for all
    if topic is None:
        topic = [False] * len(modalities)
    # check if modalities and topic lists have the same length
    if len(modalities) != len(topic):
        raise ValueError(f"Modalities and topic lists must have the same length.")
    # check if modal_weights are provided, otherwise use equal weights
    if modal_weights is None:
        modal_weights = [1] * len(modalities)
    # check if modal_weights and modalities lists have the same length
    if len(modal_weights) != len(modalities):
        raise ValueError(f"Modalities and modal_weights lists must have the same length.")

    # Find common barcodes in provided modalities
    if column_key is None:
        barcodes = set.intersection(*(set(mdata[mod].obs.index) for mod in modalities))
    else:
        barcodes = set.intersection(*(set(mod.obs[column_key]) for mod in modalities))
    barcodes = list(barcodes)

    graphs = []
    umaps = []
    for mod, top in zip(modalities, topic):
        adata = mdata.mod[mod][barcodes]

        if top == "LSA":
            dat = TOPICMODEL(
                adata,
                n_topics=n_topics,
                binarize=False,
                smart_code="lfu",
            )
            dat.runLSA()
            cell_topic = dat.get_cell_topic()
            # get graph and umap
            umap, graph = cluster_topic(cell_topic, modularityAlg="leiden", resolution=1, nk=nK)
            adata.obsm["X_pca"] = cell_topic
            adata.obsm["X_umap"] = umap
        elif top == "LDA":
            dat = TOPICMODEL(
                adata,
                n_topics=n_topics,
                binarize=False,
                n_passes=2,
                n_workers=4,
            )
            dat.runLDA()
            cell_topic = dat.get_cell_topic()
            # get graph and umap
            umap, graph = cluster_topic(cell_topic, modularityAlg="leiden", resolution=1, nk=nK)
            adata.obsm["X_pca"] = cell_topic
            adata.obsm["X_umap"] = umap
        else:
            sc.pp.neighbors(adata, n_neighbors=nK)
            sc.tl.leiden(adata)
            sc.tl.umap(adata, min_dist=0.5, spread=5, init_pos="random", random_state=42)
            # get graph and umap
            graph = get_igraph_from_adjacency(adata.obsp["connectivities"], directed=True)
            umap = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"])
            umap.index = adata.obs.index
            umap["cluster_mod"] = adata.obs.leiden

        graphs.append(graph)
        umaps.append(umap)

    # leiden multi-layer clustering
    optimiser = la.Optimiser()
    parts = []
    for graph in graphs:
        part = la.ModularityVertexPartition(graph)
        parts.append(part)

    optimiser.optimise_partition_multiplex(parts, layer_weights=modal_weights, n_iterations=-1)
    print("Detected clusters: ", set(parts[-1].membership))
    umap[-1]["cluster_multi"] = parts[-1].membership

    # merge modality UMAPs
    multi_umap = umaps[0].copy()

    for i, umap in enumerate(umaps[1:], start=1):
        multi_umap = multi_umap.merge(umap, left_index=False, right_index=True, left_on=column_key)
    multi_umap["cluster_mod"] = [int(x) for x in multi_umap["cluster_mod"].to_list()]

    return multi_umap, mdata


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
