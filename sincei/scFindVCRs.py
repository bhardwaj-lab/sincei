#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse

import numpy as np
import pandas as pd
import anndata as ad
import ruptures as rpt

from sincei import ParserCommon
from sincei.VCRfinder import VCRfinder


def parseArguments():
    io_args = ParserCommon.inputOutputOptions(opts=["h5adfile"], requiredOpts=["h5adfile"])
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        parents=[io_args, get_args(), other_args],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""
        ``scFindVCRs`` calls variable chromatin regions (VCRs) from binned chromatin data. It takes a
        .h5ad file containing single-cell genomic signal in bins, and outputs BED files with genome
        segmentations for different sensitivities. This method is decsribed in [Sancho-Gómez et al. TBP].

        First, a bin-to-bin correlation matrix is computed for each chromosome.

        Then, the correlation matrix is turned into a score map by convolving a number of square
        Gaussian kernels along its main diagonal. Each kernel has a sigma calculated using. Each kernel produces a 1-D score for
        each bin, which are stacked into a matrix where each row corresponds to a kernel scale and each column to a bin.

        Finally, the PELT change-point detection algorithm is applied to the score map to identify
        regions with distinct correlation patterns. This step depends on a penalty parameter that
        controls the number of detected regions.
        """,
        usage="scFindVCRs -i binned_signal.h5ad -bs 2000 -mr 100000 -nk 20 -p 5 10 20 -o detected_VCRs.bed",
        add_help=False,
    )

    return parser


def get_args():
    parser = argparse.ArgumentParser(add_help=False)

    general = parser.add_argument_group("VCR detection options")

    general.add_argument(
        "--binSize",
        "-bs",
        type=int,
        help="The size of the bins in the input Anndata object.",
        required=True,
    )

    general.add_argument(
        "--maxRegionSize",
        "-mr",
        type=int,
        help="""
        The maximum region size to be considered, in base pairs. Larger regions may increase compute time.
        Defaults to 100 times the bin size.
        """,
        default=None,
    )

    general.add_argument(
        "--nKernels",
        "-nk",
        type=int,
        help="""
        The number of kernels to use for the score map. More kernels generally lead to a better segmentation,
        but increase the computational cost.""",
        default=20,
    )

    general.add_argument(
        "--penalties",
        "-p",
        nargs="+",
        type=float,
        help="""
        Penalty value for change-point detection. Higher values result in fewer segments. Multiple values
        can be provided (separated by space). Each penalty value will produce a separate set of regions within
        which can be seperated from the output BED file by filtering on the "score" column.
        """,
        default=[0.5, 1, 2],
    )

    general.add_argument(
        "--outFile",
        "-o",
        type=str,
        help="""
        Name of the output file (BED format) with genome segmentation result. The penalty threshold that defines
        each segment is saved in the "score" column of the BED file, and the BED file can be filtered based on this
        column to obtain non-overlapping segments.
        """,
        required=True,
    )

    general.add_argument(
        "--region",
        "-r",
        help="""
        Region of the genome to limit the operation to - this is useful when testing parameters to
        reduce the computing time. The format is chr:start:end, for example ``--region chr10`` or
        ``--region chr10:456700:891000``.
        """,
        metavar="CHR:START:END",
        required=False,
        type=ParserCommon.genomicRegion,
    )

    return parser


def main(args=None):
    args = parseArguments().parse_args(args)

    if args.maxRegionSize is None:
        args.maxRegionSize = args.binSize * 100

    adata = ad.read_h5ad(args.input)

    pen_bed_df = VCRfinder(
        adata=adata,
        binsize=args.binSize,
        max_region=args.maxRegionSize,
        n_kernels=args.nKernels,
        penalties=args.penalties,
        region=args.region,
    )

    pen_bed_df.to_csv(args.outFile, sep="\t", header=False, index=False)

#    for pen in args.penalties:
#        out_bed_df = pen_bed_df[pen_bed_df["penalty"] == pen][["chrom", "start", "end"]]
#        out_bed_df.to_csv(f"{args.outFile}_pen{pen}.bed", sep="\t", header=False, index=False)
