#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import pandas as pd
from scipy import sparse, io
from itertools import compress

# plotting
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

# clustering and umap
from scipy.sparse import issparse, coo_matrix, csr_matrix
from sklearn.metrics import pairwise_distances
# topic models
from gensim import corpora, matutils, models
# Louvain clustering and UMAP
from networkx import convert_matrix
import leidenalg as la
import community
import umap

# single-cell stuff

import anndata
import scanpy as scp
from scanpy.neighbors import _compute_connectivities_umap,  _get_indices_distances_from_dense_matrix
from scanpy._utils import get_igraph_from_adjacency

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
                          metavar='LOOM',
                          help='Input file in the loom format',
                          required=True)

    required.add_argument('--outFile', '-o',
                         type=str,
                         required=True,
                         help='The file to write results to. The output file is an updated .loom object containing cell metadata, UMAP coordinates and cluster IDs.')

    general = parser.add_argument_group('General arguments')

    general.add_argument('--outFileUMAP', '-op',
                         type=str,
                         required=False,
                         help='The output plot file (for UMAP). If you specify this option, another file with the same '
                         'prefix (and .txt extention) is also created with the raw UMAP coordinates.')

    general.add_argument('--outFileTrainedModel', '-om',
                         type=argparse.FileType('w'),
                         required=False,
                         help='The output file for the trained LSI model. The saved model can be used later to embed/compare new cells '
                              'to the existing cluster of cells.')

    general.add_argument('--outGraph', '-og',
                         type=argparse.FileType('w'),
                         required=False,
                         help='The output file for the Graph object (lgl format) which can be used for further clustering/integration.')

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

    general.add_argument('--whitelist', '-w',
                         default=None,
                         type=argparse.FileType('r'),
                         help='A list of barcodes to be included for the clustering. The barcodes '
                         '(along with sample labels) must be present in the anndata object (Default: %(default)s)')

    general.add_argument('--method', '-m',
                         type=str,
                         choices=['LSA'],
                         default='LSA',
                         help='The dimentionality reduction method for clustering. (Default: %(default)s)')

    general.add_argument('--nPrinComps', '-n',
                         default=20,
                         type=int,
                         help='Number of principle components to reduce the dimentionality to. '
                              'Use higher number for samples with more expected heterogenity. (Default: %(default)s)')

    general.add_argument('--nNeighbors', '-nk',
                         default=30,
                         type=int,
                         help='Number of nearest neighbours to consider for clustering and UMAP. This number should be chosen considering '
                              'the total number of cells and expected number of clusters. Smaller number will lead to more fragmented clusters. '
                              '(Default: %(default)s)')

    general.add_argument('--clusterResolution', '-cr',
                         default=1.0,
                         type=float,
                         help='Resolution parameter for clustering. Values lower than 1.0 would result in less clusters, '
                              'while higher values lead to splitting of clusters. In most cases, the optimum value would be between '
                              '0.8 and 1.2. (Default: %(default)s)')

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


def main(args=None):
    args = parseArguments().parse_args(args)

    adata = anndata.read_loom(args.input)
    adata = preprocess_adata(adata, args.minCellSum, args.minRegionSum)
    #adat = lsa_anndata(adat, args.nPrinComps, args.scaleFactor)
    #adat = UMAP_clustering(adat)

    ## LSA and clustering based on gensim
    mtx = sparse.csr_matrix(adata.X.transpose())
    corpus_lsi, cell_topic, corpus_tfidf = LSA_gensim(mtx, list(adata.obs.index), list(adata.var.index), nTopics = args.nPrinComps, smartCode='lfu')
    #umap_lsi, graph = cluster_LSA(cell_topic, modularityAlg='leiden', resolution=args.clusterResolution, nk=args.nNeighbors)

    ## update the anndata object, drop cells which are not in the anndata
    #adata=adata[umap_lsi.index]
    adata.obsm['X_pca']=np.asarray(cell_topic.iloc[:,1:args.nPrinComps])
    #adata.obsm['X_umap']=np.asarray(umap_lsi.iloc[:,0:2])
    #adata.obs['cluster_lsi'] = [str(cl) for cl in umap_lsi['cluster']]
    #tfidf_mat = matutils.corpus2dense(corpus_tfidf, num_terms=len(corpus_tfidf.obj.idfs))
    #adata.layers['tfidf']=tfidf_mat.transpose()
    adata.write_loom(args.outFile, write_obsm_varm=True)

    if args.outFileUMAP:
        ## plot UMAP
        plt.rcParams['font.size'] = 8.0
        # convert cm values to inches
        fig = plt.figure(figsize=(args.plotWidth / 2.54, args.plotHeight / 2.54))
        fig.suptitle('LSA-UMAP', y=(1 - (0.06 / args.plotHeight)))
        plt.scatter(umap_lsi.UMAP1, umap_lsi.UMAP2, s=5, alpha = 0.8, c=[sns.color_palette()[x] for x in list(umap_lsi.cluster)])
        plt.tight_layout()
        plt.savefig(args.outFileUMAP, dpi=200, format=args.plotFileFormat)
        plt.close()
        prefix=args.outFileUMAP.split(".")[0]
        umap_lsi.to_csv(prefix+".tsv", sep = "\t")

    # save if asked
    if args.outFileTrainedModel:
        corpus_lsi.save(args.outFileTrainedModel)
    if args.outGraph:
        graph.write_lgl(args.outGraph)

    return 0

if __name__ == "__main__":
    main()
