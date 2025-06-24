#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse

# logs
import warnings
import logging

logger = logging.getLogger()
warnings.simplefilter(action="ignore", category=FutureWarning)

# single-cell stuff
import anndata as ad
import mudata as md
import scanpy as sc

from sincei import ParserCommon
from sincei.ParserCommon import smartLabel


def parseArguments():
    io_args = ParserCommon.inputOutputOptions(opts=["h5adfiles"])
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        parents=[io_args, get_args(), other_args],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""
        This tool combines multiple count matrices (output of scCountReads) into one, either assuming they are different samples (multi-sample)
        or different measurements on the same set of cells (multi-modal). The result is a .h5ad (AnnData) file with combined counts in
        multi-sample mode or a .h5mu (MuData) file in multi-modal mode. NOTE: it doesn't perform any 'batch effect correction' or 'integration'
        of data from different technologies, which requires more sophisticated methods.
        """,
        usage="Example usage: scCombineCounts -i sample1.h5ad sample2.h5ad -o combined.h5ad -m multi-sample > log.txt",
        add_help=False,
    )

    return parser


def get_args():
    parser = argparse.ArgumentParser(add_help=False)

    general = parser.add_argument_group("General Options")

    general.add_argument(
        "--input",
        "-i",
        metavar="H5AD",
        help="Input files in .h5ad format",
        nargs="+",
        required=True,
    )

    general.add_argument(
        "--modalities",
        "-md",
        metavar="modalities",
        help="This option is used only in multi-modal mode to specify the labels to store "
        "each of the input .h5ad files in the output .h5mu file. Multiple modalities have "
        "to be separated by a space, e.g. --modalities RNA ATAC CHiC",
        nargs="+",
        type=str,
        required=False,
    )

    general.add_argument(
        "--outFile",
        "-o",
        type=str,
        help="The file to write results to. For method: `multi-sample`, the output "
        "file is an updated .h5ad object, which can be used by other tools. "
        "For method: `multi-modal`, the output file is an .h5mu file. This file can only be "
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
    if not args.verbose:
        logger.setLevel(logging.CRITICAL)
        warnings.filterwarnings("ignore")

    if args.labels and len(args.input) != len(args.labels):
        print("The number of labels does not match the number of input files.")
        sys.exit(1)
    if not args.labels:
        # try smartlabel
        args.labels = [smartLabel(x) for x in args.input]
    adata_list = [sc.read_h5ad(x) for x in args.input]

    if args.method == "multi-sample":
        ## concatenate labels and match chrom, start, end
        var_list = []
        var_cols = ["chrom", "start", "end"]
        for lab, adata in zip(args.labels, adata_list):
            obs = adata.obs_names.to_list()
            lab = [lab] * len(obs)
            new_names = ["_".join([x, y]) for x, y in zip(lab, obs)]
            adata.obs_names = new_names
            hasinfo = all([x in adata.var.columns for x in var_cols])
            var_list.append(hasinfo)

        ## keep the chrom, start, end from original sample if present
        adata = ad.concat(adata_list)
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
        adata.write_h5ad(args.outFile)

    elif args.method == "multi-modal":
        adata_dict = dict(zip(args.modalities, adata_list))
        mdata = md.MuData(adata_dict)

        sys.stdout.write("Combined modalities: {} \n".format(len(mdata.mod)))
        sys.stdout.write("Combined cells: {} \n".format(mdata.shape[0]))
        sys.stdout.write("Combined features: {} \n".format(mdata.shape[1]))
        mdata.write_h5mu(args.outFile)

    else:
        sys.stderr.write("Choose a valid method: multi-sample or multi-modal.\n")
        sys.exit(1)

    return 0
