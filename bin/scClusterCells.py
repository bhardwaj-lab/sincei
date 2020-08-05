#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import pandas as pd
from scipy import sparse, io
from itertools import compress
from deeptools import parserCommon
# plotting
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt

# clustering and umap
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import AgglomerativeClustering
import umap

# single-cell stuff
import scanpy as scp
import anndata



def parseArguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
        This tool clusters the cells based on the input count matrix (output of scCountReads) and returns a
        tsv file with UMAP coordinates and corresponding cluster id for each barcode.
        """,
        usage='Example usage: scClusterCells.py -i cellCounts.h5ad -o clusters.tsv > log.txt')

    required = parser.add_argument_group('Required arguments')
    required.add_argument('--input', '-i',
                          metavar='ADATA',
                          help='Input count matrix, stored in h5ad format',
                          required=True)

    required.add_argument('--outFile', '-o',
                         type=parserCommon.writableFile,
                         required=True,
                         help='The file to write results to. The output file contains cell metadata, UMAP coordinates and cluster IDs.')

    general = parser.add_argument_group('General arguments')

    general.add_argument('--plotFile', '-p',
                         type=parserCommon.writableFile,
                         required=False,
                         help='The output plot file (for UMAP)')

    general.add_argument('--method', '-m',
                         type=str,
                         choices=['LSA'],
                         default='LSA',
                         help='The dimentionality reduction method for clustering. (Default: %(default)s)')

    general.add_argument('--minCellSum', '-c',
                         default=1000,
                         type=float,
                         help='For filtering of cells: minimum number of regions detected in a cell for '
                               'the cell to be kept. (Default: %(default)s)')

    general.add_argument('--minRegionSum', '-r',
                         default=100,
                         type=float,
                         help='For filtering of regions: Minimum number of cells the regions should be present in, '
                              'for the region to be kept. (Default: %(default)s)')

    general.add_argument('--scaleFactor', '-s',
                         default=100000,
                         type=float,
                         help='The scale factor to multiply with, for calculating term frequency (when --method LSA). '
                              '(Default: %(default)s)')

    general.add_argument('--nPrinComps', '-n',
                         default=20,
                         type=int,
                         help='Number of principle components to reduce the dimentionality to. '
                              'Use higher number for samples with more heterogenity. (Default: %(default)s)')

    general.add_argument('--plotWidth',
                         default=10,
                         type=float,
                         help='Output plot width (in cm). (Default: %(default)s)')

    general.add_argument('--plotHeight',
                         default=10,
                         type=float,
                         help='Output plot height (in cm). (Default: %(default)s)')

    general.add_argument('--plotFileFormat',
                          metavar='FILETYPE',
                          type=str,
                          choices=['png', 'pdf', 'svg', 'eps'],
                          default='png',
                          help='Image format type. If given, this option '
                          'overrides the image format based on the plotFile '
                          'ending. (Default: %(default)s)')

    general.add_argument('--verbose', '-v',
                         help='Set to see processing messages.',
                         action='store_true')

#    general.add_argument('--version', action='version',
#                         version='%(prog)s {}'.format(__version__))


    return parser



def read_mtx(prefix):
    mtx = io.mmread(infolder + ".counts.mtx")
    mtx = mtx.tocsr()
    with open(infolder + ".colnames.txt") as col:
        colnames = col.read().splitlines()
    col.close()
    with open(infolder + ".rownames.txt") as row:
        rownames = row.read().splitlines()
    row.close()

    return mtx, rownames, colnames

# from mtx
def preprocess_mtx(sparse_mtx, rownames, colnames, min_cell_sum, min_region_sum):
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
    scp.pp.calculate_qc_metrics(adata, inplace=True)
    # filter
    adata = adata[adata.obs.n_genes_by_counts >= min_cell_sum, :]
    adata = adata[:, adata.var.n_cells_by_counts >= min_region_sum]

    return adata

def runLSA(sparse_mtx, nPCs, scaleFactor):
    #tf = (sparse_mtx.transpose() / sparse_mtx.sum(axis=0)).transpose()
    tf = sparse_mtx / sparse_mtx.sum(axis=0)
    tf = np.log1p(tf * scaleFactor)

    idf = np.log1p(sparse_mtx.shape[1] / sparse_mtx.sum(axis=1))
    tfidf = np.multiply(tf, idf)
    svd = TruncatedSVD(n_components=nPCs)
    pca = svd.fit(np.nan_to_num(tfidf, nan = 0.0))

    return pca, tfidf

## for anndata
def lsa_anndata(adata, n_pcs, scale_factor):
    mtx = sparse.csr_matrix(adata.X)
    lsa_out, tfidf = runLSA(mtx.transpose(), n_pcs, scale_factor)

    adata.X = np.squeeze(np.asarray(tfidf.transpose()))
    adata.obsm['X_pca'] = lsa_out.components_.transpose()
    adata.uns['pca_variance'] = lsa_out.explained_variance_
    adata.uns['pca_variance_ratio'] = lsa_out.explained_variance_ratio_
    adata.layers['raw_counts'] = mtx

    return adata

def UMAP_clustering(adata):

    # make umap on PCA
    lsa_out = adata.obsm['X_pca'].transpose()[1:, :]
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, spread=1.0, metric='euclidean', init = 'random')
    embeddings = reducer.fit_transform(lsa_out.transpose())
    adata.obsm['X_umap'] = embeddings

    # louvain (from scanpy)
    scp.pp.neighbors(adata)
    scp.tl.louvain(adata)
    cluster_id = [int(x) for x in adata.obs['louvain'].to_list()]

    # return
    return adata


def main(args=None):
    args = parseArguments().parse_args(args)

    adat = anndata.read_h5ad(args.input)
    adat = preprocess_adata(adat, args.minCellSum, args.minRegionSum)
    adat = lsa_anndata(adat, args.nPrinComps, args.scaleFactor)
    adat = UMAP_clustering(adat)

    cluster_id = adat.obs.louvain.to_list()
    cluster_id = [int(x) for x in cluster_id]
    embeddings = adat.obsm['X_umap']

    if args.plotFile:
        ## plot UMAP
        plt.rcParams['font.size'] = 8.0
        # convert cm values to inches
        fig = plt.figure(figsize=(args.plotWidth / 2.54, args.plotHeight / 2.54))
        fig.suptitle('LSA-UMAP', y=(1 - (0.06 / args.plotHeight)))
        plt.scatter(
            embeddings[:, 0],
            embeddings[:, 1],
            c=[sns.color_palette()[x] for x in cluster_id])
        plt.tight_layout()
        plt.savefig(args.plotFile, dpi=200, format=args.plotFileFormat)
        plt.close()

    ## Write output clusters
    df = pd.DataFrame({'UMAP1': embeddings[:, 0],
                  'UMAP2': embeddings[:, 1],
                  'Cluster': cluster_id,
                 })

    df.index = adat.obs.index
    df.to_csv(args.outFile, sep = "\t")


    return 0

if __name__ == "__main__":
    main()
