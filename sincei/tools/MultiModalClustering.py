from __future__ import annotations

from typing import TYPE_CHECKING, cast

import leidenalg as la
import numpy as np
import pandas as pd
import scanpy as sc
import umap
from scanpy._utils import get_igraph_from_adjacency

# own modules
from sincei.tools.TopicModels import LDA, LSA
from sincei.utils import depth_components

if TYPE_CHECKING:
    from collections.abc import Sequence

    import anndata as ad
    import mudata as md
    from scipy import sparse

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

METHODS = ("PCA", "glmPCA", "LSA", "LDA")


def MultiModalClustering(
    mdata: md.MuData,
    modalities: Sequence[str] | None = None,
    method: str | Sequence[str] = "PCA",
    modal_weights: Sequence[float] | None = None,
    column_key: str | None = None,
    nK: int = 30,
    nPrinComps: int | Sequence[int] = 20,
    clusterResolution: float | Sequence[float] = 1.0,
    binarize: bool = False,
    log_normalize: bool | Sequence[bool] = False,
    drop_depth: bool | Sequence[bool] = False,
    glmPCAfamily: str = "poisson",
) -> None:
    r"""
    Performs multi-graph clustering on matched keys (barcodes) of a mudata object and
    stores the clustering results in mdata.obs["cluster_multi"]. It also stores the
    UMAP coordinates for each of the specified modalities in
    mdata[mod].obsm["X_umap"], where mod is the modality.

    Parameters
    ----------
    mudata : MuData
        MuData object containing several data modalities
    modalities : list[str]
        List of modalities to use for clustering, e.g. ["RNA", "ATAC", "ChIC"]
    method : list[str]
        What processing method to use for each modality. Choose between "PCA",
        "glmPCA", "LSA" or "LDA". Default is "PCA" for all modalities.
    modal_weights : list[float]
        Weights for each modality in the clustering process. Default is equal
        weighting. E.g. for RNA and ChIC, use [2, 1].
    column_key : str, optional
        Column name for the barcode. If None, the index of .obs for each modality is
        used.
    nK : int
        Number of nearest neighbours to consider for clustering and UMAP. This number
        should be chosen considering the total number of cells and expected number of
        clusters. Smaller number will lead to more fragmented clusters.
    nPrinComps : int or list[int]
        Number of principal components (for PCA or glmPCA) or number of topics (for
        LSA and LDA) to use for model. Use higher number for samples with more expected
        heterogenity. If list is provided, it must contain a value for each modality.
        Default is 20.
    clusteResolution : float or list[float]
        Resolution parameter for clustering. Values lower than 1.0 result in less
        clusters, while higher values lead to splitting of clusters. In most cases, the
        optimum value would be between 0.8 and 1.2. If list is provided, it must
        contain a value for each modality. Default is 1.0 .
    binarize : bool
        Whether to binarize the counts per region before dimensionality reduction (only
        for LSA/LDA).
    log_normalize : bool or list[bool]
        Whether to scale the counts of every cell to a total of 10,000 and then apply
        the log1p transform, before dimensionality reduction. Not valid with LDA, which
        needs counts. If list is provided, it must contain a value for each modality.
        Default is False.
    drop_depth : bool or list[bool]
        Whether to fit one extra component, and drop every component whose scores
        correlate with the total counts of the cells (|Spearman rho| above
        DEPTH_CORRELATION_CUTOFF). If list is provided, it must contain a value for
        each modality. Default is False.
    glmPCAfamily : str
        The choice of exponential family distribution to use for glmPCA method. Default
        is "poisson".
    """

    # check if modalities are provided, otherwise use all
    if modalities is None:
        msg = "Choose modalities to use for clustering."
        raise ValueError(msg)
    # check if modalities are in mudata object
    for mod in modalities:
        if mod not in mdata.mod:
            msg = f"Modality {mod} not found in MuData object."
            raise ValueError(msg)
    # check if method is a string, then convert to list
    methods = [method] * len(modalities) if isinstance(method, str) else list(method)
    # check if modalities and method lists have the same length
    if len(modalities) != len(methods):
        msg = "Modalities and method lists must have the same length."
        raise ValueError(msg)
    # check if modal_weights are provided, otherwise use equal weights
    if modal_weights is None:
        modal_weights = [1] * len(modalities)
    # check if modal_weights and modalities lists have the same length
    if len(modal_weights) != len(modalities):
        msg = "Modalities and modal_weights lists must have the same length."
        raise ValueError(msg)
    # check if nPrinComps are provided, otherwise use default
    if np.ndim(nPrinComps) == 0:
        n_comps = [cast("int", nPrinComps)] * len(modalities)
    else:
        n_comps = list(cast("Sequence[int]", nPrinComps))
    # check if nPrinComps and modalities lists have the same length
    if len(n_comps) != len(modalities):
        msg = "Modalities and nPrinComps lists must have the same length."
        raise ValueError(msg)
    # a single resolution applies to every modality
    if np.ndim(clusterResolution) == 0:
        resolutions = [cast("float", clusterResolution)] * len(modalities)
    else:
        resolutions = list(cast("Sequence[float]", clusterResolution))
    if len(resolutions) != len(modalities):
        msg = "Modalities and clusterResolution lists must have the same length."
        raise ValueError(msg)
    # a single bool applies to every modality
    if np.ndim(log_normalize) == 0:
        log_normalizes = [cast("bool", log_normalize)] * len(modalities)
    else:
        log_normalizes = list(cast("Sequence[bool]", log_normalize))
    if len(log_normalizes) != len(modalities):
        msg = "Modalities and log_normalize lists must have the same length."
        raise ValueError(msg)
    if np.ndim(drop_depth) == 0:
        drop_depths = [cast("bool", drop_depth)] * len(modalities)
    else:
        drop_depths = list(cast("Sequence[bool]", drop_depth))
    if len(drop_depths) != len(modalities):
        msg = "Modalities and drop_depth lists must have the same length."
        raise ValueError(msg)
    for met, normalize in zip(methods, log_normalizes, strict=True):
        if met not in METHODS:
            msg = f"{met!r} is not a valid method. Use one of {', '.join(METHODS)}."
            raise ValueError(msg)
        if met == "LDA" and normalize:
            msg = "log_normalize is not valid with LDA, which needs counts."
            raise ValueError(msg)

    # Find common barcodes in provided modalities
    if column_key is None:
        common = set.intersection(
            *(set(mdata.mod[mod].obs.index) for mod in modalities)
        )
    else:
        common = set.intersection(
            *(set(mdata.mod[mod].obs[column_key]) for mod in modalities)
        )
    barcodes = list(common)

    adatas = []
    graphs = []
    for i, (mod, met) in enumerate(zip(modalities, methods, strict=True)):
        adata = cast("ad.AnnData", mdata.mod[mod][barcodes])

        reduction = f"X_{met.lower()}"
        if drop_depths[i]:
            n_comps[i] += 1
        matrix = adata
        if log_normalizes[i]:
            matrix = adata.copy()
            sc.pp.normalize_total(matrix, target_sum=1e4)
            sc.pp.log1p(matrix)

        if met == "PCA":
            adata.obsm[reduction] = sc.pp.pca(matrix.X, n_comps=n_comps[i])

        elif met == "LSA":
            lsa = LSA(n_pc=n_comps[i], binarize=binarize)
            adata.obsm[reduction] = lsa.fit_transform(matrix)

        elif met == "LDA":
            lda = LDA(n_topics=n_comps[i], binarize=binarize)
            adata.obsm[reduction] = lda.fit_transform(matrix)

        elif met == "glmPCA":
            # import glmPCA (not imported on top due to special optional import of
            # mctorch)
            from sincei.tools.GLMPCA import GLMPCA  # noqa: PLC0415

            ## glmPCA using mctorch
            glmpca = GLMPCA(
                n_pc=n_comps[i],
                family=glmPCAfamily,
            )
            glmpca.fit(matrix)
            assert glmpca.saturated_loadings_ is not None
            cell_pcs = glmpca.saturated_loadings_.detach().numpy()

            ## update the anndata object
            adata.obsm[reduction] = np.asarray(cell_pcs)

        if drop_depths[i]:
            depth = np.asarray(adata.X.sum(axis=1), dtype=np.float64).ravel()
            keep = depth_components(adata.obsm[reduction], depth)
            if not keep.any():
                msg = f"drop_depth dropped every component of modality {mod}."
                raise ValueError(msg)
            adata.obsm[reduction] = adata.obsm[reduction][:, keep]

        sc.pp.neighbors(adata, use_rep=reduction, n_neighbors=nK)
        sc.tl.leiden(adata, resolution=resolutions[i])
        sc.tl.paga(adata)
        sc.pl.paga(adata, plot=False, threshold=0.1)
        sc.tl.umap(adata, min_dist=0.1, spread=5, init_pos="paga")

        # get graph
        connectivities = cast("sparse.csr_matrix", adata.obsp["connectivities"])
        graphs.append(get_igraph_from_adjacency(connectivities, directed=True))
        adatas.append(adata)

    # leiden multi-layer clustering
    optimiser = la.Optimiser()
    parts = [
        la.RBConfigurationVertexPartition(
            graph,
            resolution_parameter=resolutions[i],
        )
        for i, graph in enumerate(graphs)
    ]

    optimiser.optimise_partition_multiplex(
        parts, layer_weights=modal_weights, n_iterations=-1
    )
    print("Detected clusters: ", set(parts[0].membership))
    groups = np.array(parts[0].membership)
    mdata.obs["cluster_multi"] = pd.Categorical(values=groups.astype("U"))

    mods = cast("dict[str, ad.AnnData]", mdata.mod)
    for mod, adata in zip(modalities, adatas, strict=True):
        mods[mod] = adata


def umap_aligned(
    mdata: md.MuData,
    modalities: Sequence[str] | None = None,
    column_key: str | None = None,
    nK: int = 30,
    distance_metric: str = "euclidean",
) -> None:
    r"""
    Aligns the UMAP embeddings of the selected modalities in a mudata object using the
    UMAP AlignedUMAP class and stores them in mdata[mod].obsm["X_umap_aligned"], where
    mod is the modality. This produces an aligned UMAP for each modality, since the
    alignment for each may be slightly different.

    Parameters
    ----------
    mudata : MuData
        MuData object containing several data modalities
    modalities : list[str]
        List of modalities to use for clustering, e.g. ["RNA", "ATAC", "ChIC"]
    column_key : str, optional
        Column name for the barcode. If None, the index of .obs for each modality is
        used.
    nK : int
        Number of nearest neighbors to use for UMAP
    distance_metric : str
        Distance metric to use for UMAP, e.g. "euclidean", "cosine", etc.
    """
    # check if modalities are provided, otherwise use all
    if modalities is None:
        msg = "Choose modalities to use to align UMAP."
        raise ValueError(msg)
    # check if modalities are in mudata object
    for mod in modalities:
        if mod not in mdata.mod:
            msg = f"Modality {mod} not found in MuData object."
            raise ValueError(msg)

    # Find common barcodes in provided modalities
    if column_key is None:
        common = set.intersection(
            *(set(mdata.mod[mod].obs.index) for mod in modalities)
        )
    else:
        common = set.intersection(
            *(set(mdata.mod[mod].obs[column_key]) for mod in modalities)
        )
    barcodes = list(common)

    umaps = []
    for mod in modalities:
        adata = cast("ad.AnnData", mdata.mod[mod][barcodes])
        try:
            um = adata.obsm["X_umap"]
        except KeyError as exc:
            msg = (
                f"UMAP coordinates for modality {mod} not found. Please run UMAP first."
            )
            raise KeyError(msg) from exc

        umaps.append(um)

    # AlignedUMAP requires a mapping of relations between modalities.
    # In our case, the numerical index of each cell barcode to itself are the relations.
    relation_dict = {i: i for i in range(len(barcodes))}
    relation_dicts = [relation_dict.copy() for _ in range(len(modalities) - 1)]

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
        cast("ad.AnnData", mdata[mod]).obsm["X_umap_aligned"] = (
            aligned_umap.embeddings_[i]
        )
