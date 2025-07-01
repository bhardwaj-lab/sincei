#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import re
import argparse
import pandas as pd
import scanpy as sc

# logs
import warnings
import logging

logger = logging.getLogger()
warnings.simplefilter(action="ignore", category=FutureWarning)

from sincei import ParserCommon
from sincei.Utilities import gini


### ------ Functions ------
def filter_adata(
    adata,
    filter_region_dict=None,
    filter_cell_dict=None,
    bad_chrom=None,
    bad_cells=None,
):
    # 1. regions
    if bad_chrom:
        adata = adata[:, ~adata.var.chrom.isin(bad_chrom)]
    if filter_region_dict:
        for key in filter_region_dict.keys():
            adata = adata[
                :,
                (adata.var[key] >= filter_region_dict[key][0]) & (adata.var[key] <= filter_region_dict[key][1]),
            ]

    # 2. Cells
    if bad_cells:
        adata = adata[~adata.obs.index.isin(bad_cells)]
    if filter_cell_dict:
        for key in filter_cell_dict.keys():
            adata = adata[
                (adata.obs[key] >= filter_cell_dict[key][0]) & (adata.obs[key] <= filter_cell_dict[key][1]),
                :,
            ]

    return adata


"""
def make_plots(adata, fname=None):
    # plotting
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.use('Agg')
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['svg.fonttype'] = 'none'
    import seaborn as sns# seaborn colormaps conflicts with deeptools colormaps

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

    general.add_argument('--outPlot', '-op',
                         type=str,
                         help='The output plot file. This describes the distribution of filtering metrics pre and post filtering')
    if args.outPlot:
        make_plots(adata, fname=args.outPlot)
        if cellfilter or regionfilter or badcells or badchrom:
            make_plots(adata_filt, fname=args.outPlot+".filtered.pdf")
"""


def parseArguments():
    io_args = ParserCommon.inputOutputOptions(opts=["h5adfile", "outFile"])
    plot_args = ParserCommon.plotOptions()
    other_args = ParserCommon.otherOptions()
    parser = argparse.ArgumentParser(
        parents=[io_args, get_args(), plot_args, other_args],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
        This tool performs calculates multiple quality controls metrics on the input .h5ad file (output of scCountReads)
        and (optionally) filters the input file based on filterCellArgs/filterRegionArgs. The output is either an
        updated .h5ad object (if filtering is requested) or the filtering metrics (--outMetrics) and plots (--outPlot).
        """,
        usage="Example usage: scCountQC -i cellCounts.h5ad -om qc_metrics.tsv > log.txt",
        add_help=False,
    )

    return parser


def get_args():
    parser = argparse.ArgumentParser(add_help=False)

    general = parser.add_argument_group("Filtering arguments")
    general.add_argument(
        "--describe",
        "-d",
        action="store_true",
        help="Print a list of cell and region metrics available for QC/filtering.",
    )

    general.add_argument(
        "--outMetrics",
        "-om",
        type=str,
        help="Prefix of the output file with calculated QC metrics. If given, the cell metrics are printed in "
        "<prefix>.cell.tsv and region metrics as <prefix>.region.tsv",
    )

    general.add_argument(
        "--filterCellArgs",
        "-fc",
        type=str,
        help='List of arguments to filter cells. The format is "arg_name: minvalue, maxvalue; arg_name: minvalue, maxvalue; ...." '
        "where arg_name is the QC metric for cells present in the input h5ad file. In order to view all available "
        'cell filtering metrics, run scCountFilter with "--describe". The two arguments are supplied (minvalue, maxvalue) '
        "they are used as lower and upper bounds to filter cells. Make sure they are float/integer numbers.",
    )

    general.add_argument(
        "--filterRegionArgs",
        "-fr",
        type=str,
        help='List of arguments to filter regions. The format is "arg_name: minvalue, maxvalue; arg_name: minvalue; ...." '
        "where arg_name is the QC metric for regions present in the input h5ad file. In order to view all available "
        'cell filtering metrics, run scCountFilter with "--describe". The two arguments are supplied (minvalue, maxvalue) '
        "they are used as lower and upper bounds to filter cells. Make sure they are float/integer numbers.",
    )

    general.add_argument(
        "--cell_blacklist",
        "-cb",
        default=None,
        type=argparse.FileType("r"),
        help="A list of barcodes to be included for the clustering. The barcodes "
        "(along with sample labels) must be present in the input object.",
    )

    general.add_argument(
        "--chrom_blacklist",
        "-chb",
        default=None,
        type=str,
        help="A comma separated list of chromosomes to exclude. eg. chrM, chrUn",
    )

    return parser


def main(args=None):
    args = parseArguments().parse_args(args)
    if not args.verbose:
        logger.setLevel(logging.CRITICAL)
        warnings.filterwarnings("ignore")
    try:
        adata = sc.read_h5ad(args.input)  # , obs_names="obs_names", var_names="var_names")
    except:
        sys.stderr.write("\n Error: Input file can not be read (doesn't appear to be a valid anndata object) \n")
        exit()
    ## add QC stats to the anndata object
    # 1. scanpy metrics # fraction of regions/genes with signal are included in the metrics (pct_dropouts/n_genes_by_counts)
    try:
        sc.pp.calculate_qc_metrics(adata, inplace=True)
    except IndexError:  # not enough genes/regions
        sys.stderr.write("\n Error: Too few regions in the input file to perform QC \n")
        exit()

    # 2. sincei metrics
    gini_list = [gini(i, adata.X) for i in range(adata.shape[0])]
    adata.obs["gini_coefficient"] = gini_list

    if args.outMetrics:
        mat = re.sub(".txt|.tsv|.csv", "", args.outMetrics)
        obs_tsv = mat + ".cells.tsv"
        var_tsv = mat + ".regions.tsv"
        adata.obs.to_csv(obs_tsv, sep="\t", index_label="Cell_ID")
        adata.var.to_csv(var_tsv, sep="\t", index_label="Feature_ID")

    # if --describe is asked, only print the numeric vars and obs columns
    if args.describe:
        is_num_col = [(pd.api.types.is_float_dtype(x) | pd.api.types.is_integer_dtype(x)) for x in adata.obs.dtypes]
        cols = adata.obs.loc[:, is_num_col]
        sys.stdout.write("\n Cell metrics: \n")
        sys.stdout.write("Total cells: {} \n".format(adata.shape[0]))
        print(pd.DataFrame({"min": cols.min(), "max": cols.max()}))

        is_num_row = [(pd.api.types.is_float_dtype(x) | pd.api.types.is_integer_dtype(x)) for x in adata.var.dtypes]
        rows = adata.var.loc[:, is_num_row]
        sys.stdout.write("\n Feature metrics: \n")
        sys.stdout.write("Total features: {} \n".format(adata.shape[1]))
        print(pd.DataFrame({"min": rows.min(), "max": rows.max()}))
        exit()

    if args.filterCellArgs:
        cellfilter = dict()
        for x in args.filterCellArgs.strip().split(";"):
            key = x.strip().split(":")[0]
            v = x.strip().split(":")[1]
            value = [float(x) for x in v.strip().split(",")]
            cellfilter[key] = value
    else:
        cellfilter = None
    if args.filterRegionArgs:
        regionfilter = dict()
        for x in args.filterRegionArgs.strip().split(";"):
            key = x.strip().split(":")[0]
            v = x.strip().split(":")[1]
            value = [float(x) for x in v.strip().split(",")]
            regionfilter[key] = value
    else:
        regionfilter = None

    if args.cell_blacklist:
        ## read the barcode file
        with open(args.cell_blacklist, "r") as f:
            badcells = f.read().splitlines()
        f.close()
    else:
        badcells = None

    if args.chrom_blacklist:
        badchrom = args.chrom_blacklist.strip().split(",")
    else:
        badchrom = None

    if cellfilter or regionfilter or badcells or badchrom:
        sys.stdout.write("Applying filters \n")
        adata_filt = filter_adata(
            adata,
            filter_region_dict=regionfilter,
            filter_cell_dict=cellfilter,
            bad_chrom=badchrom,
            bad_cells=badcells,
        )
        sys.stdout.write("Remaining cells: {} \n".format(adata_filt.shape[0]))
        sys.stdout.write("Remaining features: {} \n".format(adata_filt.shape[1]))
        adata_filt.write_h5ad(args.outFile)

    return 0
