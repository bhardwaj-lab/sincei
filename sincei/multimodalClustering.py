import numpy as np
import umap
import pandas as pd
import scanpy as sc

from scanpy._utils import get_igraph_from_adjacency
import leidenalg as la

# own modules
from sincei.TopicModels import TOPICMODEL
from sincei.ParserCommon import numberOfProcessors

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


def multiModal_clustering(mdata, modalities=None, method="PCA", modal_weights=None, column_key=None, nK=30,
                          nPrinComps=20, clusterResolution=1.0, binarize=False, glmPCAfamily="poisson"):
    r"""
    Performs multi-graph clustering on matched keys(barcodes) of a mudata object and stores the clustering results
    in mdata.obs["cluster_multi"]. It also stores the UMAP coordinates for each of the specified modalities in
    mdata[mod].obsm["X_umap"], where mod is the modality.
    Note: If method is "PCA" or "logPCA", the data matrix of the modality will be normalized, and log1p-transformed
    in the case of logPCA.

    Parameters
    ----------
    mudata : MuData
        MuData object containing several data modalities
    modalities : list[str]
        List of modalities to use for clustering, e.g. ["RNA", "ATAC", "CHiC"]
    method : list[str]
        What processing method to use for each modality. Choose between "PCA", "logPCA", "glmPCA", "LSA" or "LDA".
        Default is "PCA" for all modalities.
    modal_weights : list[float]
        Weights for each modality in the clustering process. Default is equal weighting. E.g. for RNA and CHiC, use [2, 1].
    column_key : str, optional
        Column name for the barcode. If None, the index of .obs for each modality is used.
    nK : int
        Number of nearest neighbours to consider for clustering and UMAP. This number should be chosen considering 
        the total number of cells and expected number of clusters. Smaller number will lead to more fragmented clusters.
    nPrinComps : int or list[int]
        Number of principal components (for logPCA or glmPCA) or number of topics (for LSA and LDA) to use for model.
        Use higher number for samples with more expected heterogenity. If list is provided, it must contain a value for each
        modality. Default is 20.
    clusteResolution : float
        Resolution parameter for clustering. Values lower than 1.0 result in less clusters, while higher values lead to
        splitting of clusters. In most cases, the optimum value would be between 0.8 and 1.2. Default is 1.0 .
    binarize : bool
        Whether to binarize the counts per region before dimensionality reduction (only for LSA/LDA).
    glmPCAfamily : str
        The choice of exponential family distribution to use for glmPCA method. Default is "poisson".
    """

    # check if modalities are provided, otherwise use all
    if modalities is None:
        raise ValueError(f"Choose modalities to use for clustering.")
    # check if modalities are in mudata object
    for mod in modalities:
        if mod not in mdata.mod.keys():
            raise ValueError(f"Modality {mod} not found in MuData object.")
    # check if method is provided, otherwise use "PCA" for all
    if method == "PCA":
        method = ["PCA"] * len(modalities)
    # check if modalities and method lists have the same length
    if len(modalities) != len(method):
        raise ValueError(f"Modalities and method lists must have the same length.")
    # check if modal_weights are provided, otherwise use equal weights
    if modal_weights is None:
        modal_weights = [1] * len(modalities)
    # check if modal_weights and modalities lists have the same length
    if len(modal_weights) != len(modalities):
        raise ValueError(f"Modalities and modal_weights lists must have the same length.")

    # Find common barcodes in provided modalities
    if column_key is None:
        barcodes = set.intersection(*(set(mdata.mod[mod].obs.index) for mod in modalities))
    else:
        barcodes = set.intersection(*(set(mdata.mod[mod].obs[column_key]) for mod in modalities))
    barcodes = list(barcodes)

    adatas = []
    graphs = []
    for mod, met in zip(modalities, method):
        adata = mdata.mod[mod][barcodes]

        if met == "PCA":
            # if no method is provided, use PCA
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.pca(adata, nPrinComps)

        elif met == "logPCA":
            ## log1p+PCA using scanpy
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.pca(adata, nPrinComps)

        elif met == "LSA":
        ## LSA using gensim
            model_object = TOPICMODEL(
                adata,
                n_topics=nPrinComps,
                binarize=binarize,
                smart_code="lfu",
            )
            model_object.runLSA()
            adata.obsm["X_pca"] = model_object.get_cell_topic()

        elif met == "LDA":
        ## LDA using gensim
            model_object = TOPICMODEL(
                adata,
                n_topics=nPrinComps,
                binarize=binarize,
                n_passes=2,
                n_workers=numberOfProcessors("max"),
            )
            model_object.runLDA()
            adata.obsm["X_pca"] = model_object.get_cell_topic()

        elif met == "glmPCA":
            # import glmPCA (not imported on top due to special optional import of mctorch)
            from sincei.GLMPCA import GLMPCA

            ## glmPCA using mctorch
            model_object = GLMPCA(
                n_pc=nPrinComps,
                family=glmPCAfamily,
            )
            model_object.fit(adata)
            cell_pcs = model_object.saturated_loadings_.detach().numpy()

            ## update the anndata object
            adata.obsm["X_pca"] = np.asarray(cell_pcs)

        sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=nK)
        sc.tl.leiden(adata, resolution=clusterResolution)
        sc.tl.paga(adata)
        sc.pl.paga(adata, plot=False, threshold=0.1)
        sc.tl.umap(adata, min_dist=0.1, spread=5, init_pos="paga")

        # get graph
        graph = get_igraph_from_adjacency(adata.obsp["connectivities"], directed=True)

        adatas.append(adata)
        graphs.append(graph)

    # leiden multi-layer clustering
    optimiser = la.Optimiser()
    parts = []
    for graph in graphs:
        part = la.ModularityVertexPartition(graph)
        parts.append(part)

    optimiser.optimise_partition_multiplex(parts, layer_weights=modal_weights, n_iterations=-1)
    print("Detected clusters: ", set(parts[0].membership))
    mdata.obs["cluster_multi"] = parts[0].membership

    for mod, adata in zip(modalities, adatas):
            mdata.mod[mod] = adata


def umap_aligned(mdata, modalities=None, column_key=None, nK=30, distance_metric="euclidean"):
    r"""
    Aligns the UMAP embeddings of the selected modalities in a mudata object using the UMAP AlignedUMAP
    class and stores them in mdata[mod].obsm["X_umap_aligned"], where mod is the modality. This produces
    an aligned UMAP for each modality, since the alignment for each may be slightly different.

    Parameters
    ----------
    mudata : MuData
        MuData object containing several data modalities
    modalities : list[str]
        List of modalities to use for clustering, e.g. ["RNA", "ATAC", "CHiC"]
    column_key : str, optional
        Column name for the barcode. If None, the index of .obs for each modality is used.
    nK : int
        Number of nearest neighbors to use for UMAP
    distance_metric : str
        Distance metric to use for UMAP, e.g. "euclidean", "cosine", etc.
    """
    # check if modalities are provided, otherwise use all
    if modalities is None:
        raise ValueError(f"Choose modalities to use to align UMAP.")
    # check if modalities are in mudata object
    for mod in modalities:
        if mod not in mdata.mod.keys():
            raise ValueError(f"Modality {mod} not found in MuData object.")

    # Find common barcodes in provided modalities
    if column_key is None:
        barcodes = set.intersection(*(set(mdata.mod[mod].obs.index) for mod in modalities))
    else:
        barcodes = set.intersection(*(set(mdata.mod[mod].obs[column_key]) for mod in modalities))
    barcodes = list(barcodes)

    adatas = []
    umaps = []
    for mod in modalities:
        adata = mdata.mod[mod][barcodes]
        try:
            um = adata.obs[:, f"UMAP1_{mod}", f"UMAP2_{mod}"]
        except KeyError:
            raise KeyError(f"UMAP coordinates for modality {mod} not found. Please run UMAP first.")

        adatas.append(adata)
        umaps.append(um)

    # AlignedUMAP requires a mapping of relations between modalities.
    # In our case, the numerical index of each cell barcode to itself are the relations.
    relation_dict = {i: i for i in range(len(barcodes))}
    relation_dicts = [relation_dict.copy() for i in range(len(modalities) - 1)]

    # UMAP
    UA = umap.AlignedUMAP(
        spread=10,
        min_dist=0.01,
        n_neighbors=nK,
        metric=distance_metric,
        init="random",
        random_state=42,
    )
    aligned_umap = UA.fit(umaps, relations=relation_dicts)

    # Update the mudata object with the aligned UMAP coordinates
    for i, mod in enumerate(modalities):
        mdata[mod].obsm["X_umap_aligned"] = aligned_umap.embeddings_[i]
    