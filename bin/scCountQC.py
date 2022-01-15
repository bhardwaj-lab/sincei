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


import scanpy as scp

### ------ Functions ------
def filter_adata(adata,
                filter_region_dict=None, filter_cell_dict=None,
                bad_chrom=None, bad_cells=None):

    # 1. regions
    if bad_chrom:
        adata=adata[:, ~adata.var.chrom.isin(bad_chrom)]
    if filter_region_dict:
        for key in filter_region_dict.keys():
            adata=adata[:, (adata.var[key] >= filter_region_dict[key][0]) &
                            (adata.var[key] <= filter_region_dict[key][1]) ]

    # 2. Cells
    if bad_cells:
        adata=adata[~adata.obs.index.isin(bad_cells)]
    if filter_cell_dict:
        for key in filter_cell_dict.keys():
            adata=adata[:, (adata.obs[key] >= filter_cell_dict[key][0]) &
                            (adata.obs[key] <= filter_cell_dict[key][1]) ]

    return adata

def make_plots(adata, fname=None):
    ## filtering of regions
    plt.figure()
    pl1=sns.violinplot(data=adata.var, x='total_counts', y='chrom')
    plt.figure()
    pl2=sns.distplot(np.log10(adata.var['total_counts']+1))
    plt.figure()
    pl3=sns.distplot(adata.var['mean_counts'])
    plt.figure()
    pl4=sns.distplot(adata.var['fraction_cells_with_signal'])
    ## cells
    plt.figure()
    pl5=sns.scatterplot(data=adata.obs, x='total_counts', y='pct_counts_in_top_100_genes')
    pl5.set(xscale="log")
    plt.figure()
    pl6=sns.distplot(adata.obs['fraction_regions_with_signal'])
    plt.figure()
    pl7=sns.scatterplot(data=adata.obs, x='total_counts', y='fraction_regions_with_signal')

    plist=[pl1, pl2, pl3, pl4, pl5, pl6, pl7]
    if fname:
        with PdfPages(fname) as pp:
            for plot in plist:
                pp.savefig(plot.figure)
    return plist

def parseArguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
        This tool clusters the cells based on the input count matrix (output of scCountReads) and returns a
        tsv file with UMAP coordinates and corresponding cluster id for each barcode.
        """,
        usage='Example usage: scCountQC.py -i cellCounts.h5ad -o clusters.tsv > log.txt')

    required = parser.add_argument_group('Required arguments')
    required.add_argument('--input', '-i',
                          metavar='LOOM',
                          help='Input file in the loom format',
                          required=True)

    general = parser.add_argument_group('General arguments')
    general.add_argument('--outPlot', '-o',
                         type=str,
                         help='The output plot file. This describes the distribution of filtering metrics pre and post filtering')

    general.add_argument('--describe', '-d',
                         action='store_true',
                         help='Print a list of cell and region metrics available for QC/filtering.')

    general.add_argument('--filterCellArgs', '-fc',
                         type=str,
                         help='List of arguments to filter cells. The format is "arg_name: minvalue, maxvalue; arg_name: minvalue, maxvalue; ...." '
                         'where arg_name is the QC metric for cells present in the input loom file. In order to view all available '
                         'cell filtering metrics, run scCountFilter with "--describe". The two arguments are supplied (minvalue, maxvalue) '
                         'they are used as lower and upper bounds to filter cells. Make sure they are float/integer numbers.')

    general.add_argument('--filterRegionArgs', '-fr',
                         type=str,
                         help='List of arguments to filter regions. The format is "arg_name: minvalue, maxvalue; arg_name: minvalue; ...." '
                         'where arg_name is the QC metric for regions present in the input loom file. In order to view all available '
                         'cell filtering metrics, run scCountFilter with "--describe". The two arguments are supplied (minvalue, maxvalue) '
                         'they are used as lower and upper bounds to filter cells. Make sure they are float/integer numbers.')

    general.add_argument('--cell_blacklist', '-cb',
                         default=None,
                         type=argparse.FileType('r'),
                         help='A list of barcodes to be included for the clustering. The barcodes '
                         '(along with sample labels) must be present in the input object.')

    general.add_argument('--chrom_blacklist', '-chb',
                         default=None,
                         type=str,
                         help='A comma separated list of chromosomes to exclude. eg. chrM, chrUn')

    general.add_argument('--outFile', '-oa',
                         type=str,
                         help='The file to write results to (if filtering is requested). The output file is an updated .loom object post-filtering.')

    general.add_argument('--verbose', '-v',
                         help='Set to see processing messages.',
                         action='store_true')

#    general.add_argument('--version', action='version',
#                         version='%(prog)s {}'.format(__version__))


    return parser


def main(args=None):
    args = parseArguments().parse_args(args)

    adata = scp.read_loom(args.input)
    # 2. Gini coefficient
    #gini_list=[]
    #for i in range(adata.shape[0]):
    #    ar=adata.X[:,i].todense()
    #    ar=np.array(ar).flatten()
    #    if len(ar[ar>0]) > 2:
    #        gini_list.append(gini(ar[ar>0]))
    #    else:
    #        gini_list.append(1.0)
    #adata.obs['gini_coefficient'] = gini_list

    # if --describe is asked, only print the numeric vars and obs columns
    if args.describe:
        cols=adata.obs.loc[:,adata.obs.dtypes.isin(['int', 'float64'])]
        print("Cell metrics:")
        print(pd.DataFrame({'min': cols.min(),
                      'max': cols.max()}))

        rows=adata.var.loc[:,adata.var.dtypes.isin(['int', 'float64'])]
        print("Region metrics:")
        print(pd.DataFrame({'min': rows.min(),
                            'max': rows.max()}))
        exit()

    if args.filterCellArgs:
        cellfilter=dict()
        for x in args.filterCellArgs.strip().split(";"):
            key=x.strip().split(":")[0]
            v=x.strip().split(":")[1]
            value=[float(x) for x in v.strip().split(",")]
            cellfilter[key] = value
    else:
        cellfilter=None
    if args.filterRegionArgs:
        regionfilter=dict()
        for x in args.filterRegionArgs.strip().split(";"):
            key=x.strip().split(":")[0]
            v=x.strip().split(":")[1]
            value=[float(x) for x in v.strip().split(",")]
            regionfilter[key] = value
    else:
        regionfilter=None

    if args.cell_blacklist:
        ## read the barcode file
        with open(args.cell_blacklist, 'r') as f:
            badcells = f.read().splitlines()
        f.close()
    else:
        badcells=None

    if args.chrom_blacklist:
        badchrom=args.chrom_blacklist.strip().split(",")
    else:
        badchrom=None

    if cellfilter or regionfilter or badcells or badchrom:
        print(cellfilter)
        print(regionfilter)
        adata_filt=filter_adata(adata,
                        filter_region_dict=regionfilter,
                        filter_cell_dict=cellfilter,
                        bad_chrom=badchrom,
                        bad_cells=badcells)
        print(adata_filt.shape)
        adata_filt.write_loom(args.outFile)

    if args.outPlot:
        make_plots(adata, fname=args.outPlot)
        if cellfilter or regionfilter or badcells or badchrom:
            make_plots(adata_filt, fname=args.outPlot+".filtered.pdf")

    return 0

if __name__ == "__main__":
    main()
