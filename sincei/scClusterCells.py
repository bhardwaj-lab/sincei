#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import copy
import argparse
import numpy as np
import pandas as pd
import torch
from scipy import sparse, io
from sklearn.preprocessing import binarize

# logs
import warnings
import logging

logger = logging.getLogger()
warnings.simplefilter(action="ignore", category=FutureWarning)

# plotting
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"

# single-cell stuff
import anndata
import scanpy as sc


## own Functions
from sincei import ParserCommon
from sincei.TopicModels import TOPICMODEL

from sincei.GLMPCA import EXPONENTIAL_FAMILY_DICT  # , GLMPCA


def parseArguments():
    io_args = ParserCommon.inputOutputOptions(opts=["loomfile", "outFile"], requiredOpts=["outFile"])
    plot_args = ParserCommon.plotOptions()
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        parents=[io_args, get_args(), plot_args, other_args],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,  # argparse.RawDescriptionHelpFormatter,
        description="""
        This tool clusters the cells based on the input count matrix (output of scCountReads) and performs dimentionality reduction,
        community detection and 2D projection (UMAP) of the cells. The result is an updated loom object, and (optionally) a plot file
        and a tsv file with UMAP coordinates and corresponding cluster id for each barcode.
        """,
        usage="Example usage: scClusterCells -i cellCounts.loom -o clustered.loom -op <umap_prefix>.png  > log.txt",
        add_help=False,
    )

    return parser


def get_args():
    parser = argparse.ArgumentParser(add_help=False)

    general = parser.add_argument_group("Clustering Options")
    general.add_argument(
        "--outFileUMAP",
        "-op",
        type=str,
        required=False,
        help="The output plot file (for UMAP). If you specify this option, another file with the same "
        "prefix (and .txt extention) is also created with the raw UMAP coordinates.",
    )

    general.add_argument(
        "--outFileTrainedModel",
        "-om",
        type=argparse.FileType("w"),
        required=False,
        help="The output file for the trained LSI model. The saved model can be used later to embed/compare new cells "
        "to the existing cluster of cells.",
    )

    general.add_argument(
        "--method",
        "-m",
        type=str,
        choices=["logPCA", "LSA", "LDA", "glmPCA"],
        default="LSA",
        help="The dimentionality reduction method for clustering. (Default: %(default)s)",
    )

    general.add_argument(
        "--glmPCAfamily",
        "-gf",
        type=str,
        choices=EXPONENTIAL_FAMILY_DICT.keys(),
        default="poisson",
        help="The choice of exponential family distribution to use for glmPCA method. (Default: %(default)s)",
    )

    general.add_argument(
        "--binarize",
        action="store_true",
        help="Binarize the counts per region before dimentionality reduction (only for LSA/LDA)",
    )

    general.add_argument(
        "--nPrinComps",
        "-n",
        default=20,
        type=int,
        help="Number of principle components to reduce the dimentionality to. "
        "Use higher number for samples with more expected heterogenity. (Default: %(default)s)",
    )

    general.add_argument(
        "--nNeighbors",
        "-nk",
        default=30,
        type=int,
        help="Number of nearest neighbours to consider for clustering and UMAP. This number should be chosen considering "
        "the total number of cells and expected number of clusters. Smaller number will lead to more fragmented clusters. "
        "(Default: %(default)s)",
    )

    general.add_argument(
        "--clusterResolution",
        "-cr",
        default=1.0,
        type=float,
        help="Resolution parameter for clustering. Values lower than 1.0 would result in less clusters, "
        "while higher values lead to splitting of clusters. In most cases, the optimum value would be between "
        "0.8 and 1.2. (Default: %(default)s)",
    )

    return parser


def main(args=None):
    args = parseArguments().parse_args(args)
    if not args.verbose:
        logger.setLevel(logging.CRITICAL)
        warnings.filterwarnings("ignore")

    adata = sc.read_loom(args.input, obs_names="obs_names", var_names="var_names")
    mtx = sparse.csr_matrix(adata.X.copy().transpose())  # features x cells
    cells = copy.deepcopy(adata.obs_names.to_list())
    regions = copy.deepcopy(adata.var_names.to_list())

    if args.binarize:
        mtx = binarize(mtx, copy=True)

    if args.method == "logPCA":
        ## log1p+PCA using scanpy
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.pca(adata, args.nPrinComps)

    elif args.method == "LSA":
        ## LSA using gensim
        model_object = TOPICMODEL(
            mtx,
            cells,
            regions,
            n_topics=args.nPrinComps,
            smart_code="lfu",
        )
        model_object.runLSA()
        cell_topic = model_object.get_cell_topic(pop_sparse_cells=True)
        ## update the anndata object, drop cells which are not in the anndata, drop 1st PC
        adata = adata[cell_topic.index]
        adata.obsm["X_pca"] = np.asarray(cell_topic.iloc[:, 1 : args.nPrinComps])

    elif args.method == "LDA":
        ## LDA using gensim
        model_object = TOPICMODEL(
            mtx,
            cells,
            regions,
            n_topics=args.nPrinComps,
            n_passes=2,
            n_workers=4,
        )
        model_object.runLDA()
        cell_topic = model_object.get_cell_topic(pop_sparse_cells=True)
        ## update the anndata object, drop cells which are not in the anndata, drop 1st PC
        adata = adata[cell_topic.index]
        adata.obsm["X_pca"] = np.asarray(cell_topic.iloc[:, 1 : args.nPrinComps])

    elif args.method == "glmPCA":
        # import glmPCA (not imported on top due to special optional import of mctorch)
        from sincei.GLMPCA import GLMPCA

        # convert mtx to torch tensor
        mtx = torch.tensor(mtx.todense())  # feature*cell tensor

        ## glmPCA using mctorch
        model_object = GLMPCA(
            n_pc=args.nPrinComps,
            family=args.glmPCAfamily,
        )
        model_object.fit(mtx)
        cell_pcs = model_object.saturated_loadings_.detach().numpy()

        ## update the anndata object
        adata.obsm["X_pca"] = np.asarray(cell_pcs)

    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=args.nNeighbors)
    sc.tl.leiden(adata, resolution=args.clusterResolution)
    sc.tl.paga(adata)
    sc.pl.paga(adata, plot=False, threshold=0.1)
    sc.tl.umap(adata, min_dist=0.1, spread=5, init_pos="paga")

    adata.write_loom(args.outFile, write_obsm_varm=True)

    if args.outFileUMAP:
        ## plot UMAP
        fig = sc.pl.umap(adata, color="leiden", legend_loc="on data", return_fig=True, show=False)
        fig.savefig(args.outFileUMAP, dpi=300, format=args.plotFileFormat)
        prefix = args.outFileUMAP.split(".")[0]
        umap_lsi = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs.index)
        umap_lsi["cluster"] = adata.obs["leiden"]
        umap_lsi.to_csv(prefix + ".tsv", sep="\t", index_label="Cell_ID")
        # plt.rcParams['font.size'] = 8.0
        # convert cm values to inches
        # fig = plt.figure(figsize=(args.plotWidth / 2.54, args.plotHeight / 2.54))
        # fig.suptitle('LSA-UMAP', y=(1 - (0.06 / args.plotHeight)))
        # plt.scatter(umap_lsi.UMAP1, umap_lsi.UMAP2, s=5, alpha = 0.8, c=[sns.color_palette()[x] for x in list(umap_lsi.cluster)])
        # plt.tight_layout()
        # plt.savefig(args.outFileUMAP, dpi=200, format=args.plotFileFormat)
        # plt.close()

    # save if asked
    if args.outFileTrainedModel:
        model_object.lsi_model.save(args.outFileTrainedModel)
    #    if args.outGraph:
    # 'The output file for the Graph object (lgl format) which can be used for further clustering/integration.'
    #        graph.write_lgl(args.outGraph)

    return 0
