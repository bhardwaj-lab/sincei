#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import numpy as np
import pandas as pd
from scipy import sparse, io

# single-cell stuff
import anndata
import scanpy as sc


## own Functions
# scriptdir=os.path.abspath(os.path.join(__file__, "../../sincei"))
# sys.path.append(scriptdir)
from sincei import ParserCommon
from sincei.Clustering import LSA_gensim


def parseArguments():
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        parents=[get_args(), other_args],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""
        This tool combines multiple count matrices (output of scCountReads) into one, either assuming they are different samples (multi-sample)
        or different measurements on the same set of cells (multi-modal). The result is a .loom file with combined counts. NOTE: it doesn't perform
        any 'batch effect correction' or 'integration' of data from different technologies, which requires more sophisticated methods.
        """,
        usage="Example usage: scCombineCounts -i sample1.loom sample2.loom -o combined.loom  > log.txt",
        add_help=False,
    )

    return parser


def get_args():
    parser = argparse.ArgumentParser(add_help=False)

    general = parser.add_argument_group("General Options")

    general.add_argument(
        "--input",
        "-i",
        metavar="LOOM",
        help="Input files in .loom format",
        nargs="+",
        required=True,
    )

    general.add_argument(
        "--outFile",
        "-o",
        type=str,
        help="The file to write results to. For method: `multi-sample`, the output "
        "file is an updated .loom object, which can be used by other tools. "
        "For method: `multi-omic`, the output file is an .hdf5 file. This file can only be "
        "used by scClusterCells, to perform multi-modal clustering. ",
        required=True,
    )

    general.add_argument(
        "--labels",
        "-l",
        metavar="sample1 sample2",
        help="User defined labels instead of default labels from "
        "file names. Multiple labels have to be separated by a space, e.g. "
        "--labels sample1 sample2 sample3",
        nargs="+",
    )

    general.add_argument(
        "--method",
        "-m",
        type=str,
        choices=["multi-sample", "multi-modal"],
        default="multi-sample",
        help="How to merge the counts from the provided samples. "
        "`multi-sample`: assumes that each sample is the independent, "
        "but were counted in the same manner (i.e. on same features), therefore "
        "it looks for feature overlaps, but not for barcode overlaps. "
        "`multi-modal`: assumes that the counts were generated in 2 different ways, "
        "but from the same set of cells (for example, using a multi-omic technology), "
        "therefore it looks for the overlap of cell barcodes, but not for the overlaps "
        "of features (Default: %(default)s)",
    )

    return parser


def main(args=None):
    args = parseArguments().parse_args(args)
    if args.method != "multi-sample":
        sys.stderr.write("Only multi-sample method is currently implemented")
        sys.exit(1)

    if args.labels and len(args.input) != len(args.labels):
        print("The number of labels does not match the number of input files.")
        sys.exit(1)
    if not args.labels:
        if args.smartLabels:
            args.labels = smartLabels(args.input)
        else:
            args.labels = [os.path.basename(x) for x in args.input]
    adata_list = [
        sc.read_loom(x, obs_names="obs_names", var_names="var_names")
        for x in args.input
    ]

    ## concatenate labels and match chrom, start, end
    var_list = []
    var_cols = ["chrom", "start", "end"]
    for lab, ad in zip(args.labels, adata_list):
        obs = ad.obs_names.to_list()
        lab = [lab] * len(obs)
        new_names = ["_".join([x, y]) for x, y in zip(lab, obs)]
        ad.obs_names = new_names
        hasinfo = all([x in ad.var.columns for x in var_cols])
        var_list.append(hasinfo)

    ## keep the chrom, start, end from original sample if present
    adata = anndata.concat(adata_list)
    if all(var_list):
        var_df = adata_list[0].var[var_cols]
        adata.var = adata.var.join(var_df)
    else:
        sys.stderr.write(
            "WARNING: Not all input files contain the 'chrom', 'start', 'end' information. "
            "The output will lack these fields. This might cause an error in some downstream tools"
        )

    sys.stdout.write("Combined cells: {} \n".format(adata.shape[0]))
    sys.stdout.write("Combined features: {} \n".format(adata.shape[1]))
    adata.write_loom(args.outFile)
    return 0
