#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import os

from deeptools import parserCommon, bamHandler, utilities
from deeptools.mapReduce import mapReduce
from deeptoolsintervals import GTF
import numpy as np
import pandas as pd

# plotting
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"

# logs
import warnings
import logging

logger = logging.getLogger()
warnings.simplefilter(action="ignore", category=FutureWarning)

from sincei.Utilities import *
from sincei import ParserCommon


def parseArguments():
    io_args = ParserCommon.inputOutputOptions(opts=["bamfile", "whitelist", "outFile"], requiredOpts=["bamfile"])
    bam_args = ParserCommon.bamOptions(
        suppress_args=["labels", "smartLabels", "distanceBetweenBins", "region"],
        default_opts={"binSize": 100000},
    )
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        parents=[io_args, get_args(), bam_args, other_args],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
        This tool identifies barcodes present in a BAM files and produces a list. You can optionally filter these
        barcodes by matching them to a whitelist or based on total counts.
        """,
        usage="Example usage: scFilterBarcodes -b sample.bam -w whitelist.txt > barcodes_detected.txt",
        add_help=False,
    )

    return parser


def get_args():
    parser = argparse.ArgumentParser(add_help=False)
    general = parser.add_argument_group("Counting Options")

    general.add_argument(
        "--minHammingDist",
        "-d",
        help="Minimum hamming distance to match the barcode in whitelist. Note that increasing the "
        "hamming distance really slows down the barcode detection process.",
        metavar="INT",
        type=int,
        default=0,
        required=False,
    )

    general.add_argument(
        "--minCount",
        "-mc",
        help="Minimum no. of bins with non-zero counts, in order to report a barcode. Note that this number would range "
        "from 0 to genome size/binSize. ",
        metavar="INT",
        type=int,
        default=0,
        required=False,
    )

    general.add_argument(
        "--minMappingQuality",
        "-mq",
        metavar="INT",
        help="If set, only reads that have a mapping " "quality score of at least this are " "considered.",
        type=int,
    )

    general.add_argument(
        "--rankPlot",
        "-rp",
        type=parserCommon.writableFile,
        help='The output file name to plot the ranked counts per barcode (similar to the "knee plot",'
        "but counts in this case would be the number of non-zero bins)",
    )

    return parser


##
def ham_dist(s1, s2):
    r"""get hamming dist between 2 sequences

    Parameters
    ----------
    s1, s2 : string, sequences

    Returns
    ----------
    integer hamming distance
    """
    if len(s1) != len(s2):
        raise ValueError("Undefined")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def getFiltered_worker(arglist):
    r"""get a set of filtered barcodes based on the provided criteria

    Parameters
    ----------
    arglist : list of filtering arguments

    Returns
    ----------
    BCset: a set of barcodes
    """
    chrom, start, end, args = arglist
    # Fix the bounds
    if end <= start:
        end = start + 1
    ## open blacklist file
    blackList = None
    if args.blackListFileName is not None:
        blackList = GTF(args.blackListFileName)

    fh = bamHandler.openBam(args.bamfile)
    chromUse = utilities.mungeChromosome(chrom, fh.references)

    BCset = set()
    for read in fh.fetch(chromUse, start, end):
        if read.pos < start:
            # ensure that we never double count (in case distanceBetweenBins == 0)
            continue

        if read.flag & 4:
            # Ignore unmapped reads, they were counted already
            continue

        if args.minMappingQuality and read.mapq < args.minMappingQuality:
            continue

        ## reads in blacklisted regions
        if blackList and blackList.findOverlaps(
            chrom,
            read.reference_start,
            read.reference_start + read.infer_query_length(always=False) - 1,
        ):
            continue

        ## get barcode and append to set
        try:
            bc = read.get_tag(args.cellTag)
        except KeyError:
            continue
        if args.whitelist:
            # match barcode to whitelist
            if args.minHammingDist == 0:
                if bc in args.whitelist:
                    BCset.add(bc)
            else:
                try:
                    hamdist = [ham_dist(x, bc) for x in args.whitelist]
                    if min(hamdist) <= args.minHammingDist:
                        BCset.add(bc)
                except ValueError:
                    continue
        else:
            BCset.add(bc)
    fh.close()

    # return the set with barcodes detected
    return BCset


#
def count_occurrences(res):
    r"""count occurances of elements (barcodes) in a list of sets

    Parameters
    ----------
    res : list of sets of barcodes

    Returns
    ----------
    counts
    """
    barcodes = set.union(*res)
    counts = {barcode: 0 for barcode in barcodes}
    for set_ in res:
        for barcode in set_:
            counts[barcode] += 1
    return counts


def main(args=None):
    args = parseArguments().parse_args(args)
    if not args.verbose:
        logger.setLevel(logging.CRITICAL)
        warnings.filterwarnings("ignore")

    # open barcode file and read the content in a list
    if args.whitelist:
        with open(args.whitelist, "r") as f:
            barcodes = f.read().splitlines()
        args.whitelist = barcodes

    bhs = bamHandler.openBam(args.bamfile, returnStats=True, nThreads=args.numberOfProcessors)[0]
    chrom_sizes = list(zip(bhs.references, bhs.lengths))
    bhs.close()

    # Get the remaining metrics
    res = mapReduce(
        [args],
        getFiltered_worker,
        chrom_sizes,
        genomeChunkLength=args.binSize,
        blackListFileName=args.blackListFileName,
        numberOfProcessors=args.numberOfProcessors,
        verbose=args.verbose,
    )
    ## res, should be a list of sets
    df = pd.DataFrame.from_dict(count_occurrences(res), orient="index", columns=["count"]).reset_index()
    df.columns = ["barcode", "count"]
    df["selected"] = True

    if args.minCount:
        if args.minCount > len(res):
            print("minCount bigger than No. of bins. Reducing to maximum")
            args.minCount = len(res)
        final_set = df.loc[df["count"] >= args.minCount]["barcode"].to_list()
        df.loc[df["count"] < args.minCount, "selected"] = False
    else:
        final_set = df["barcode"].tolist()

    # convert count to log10
    if args.rankPlot:
        df["count_log10"] = np.log10(df["count"])
        df["count_rank"] = df["count"].rank(method="min", ascending=False)
        xrange = np.arange(0, np.round(max(df["count_rank"]), -3), np.round(int(max(df["count_rank"]) / 10), -3))
        yrange = np.arange(
            np.round(min(df["count_log10"]), 2),
            np.round(max(df["count_log10"]), 2),
            np.round(max(df["count_log10"]) / 10, 2),
        )

        fig, ax = plt.subplots()
        plt.plot(
            "count_rank",
            "count_log10",
            data=df,
            linestyle="none",
            marker="o",
            color="grey",
            markersize=0.5,
        )
        ax.set_xlabel("Barcode Rank", fontsize=12)
        ax.set_ylabel("No. of nonzero bins (log10)", fontsize=12)
        ax.set_title("Ranked counts (#bins) for detected Barcodes", fontsize=13)
        plt.xticks(xrange, fontsize=10)
        plt.yticks(yrange, fontsize=10)

        # Annotation
        if args.minCount:
            plt.axhline(np.log10(args.minCount), color="r")

        fig.tight_layout()
        plt.savefig(
            args.rankPlot,
            dpi="figure",
            format=None,
            metadata=None,
            pad_inches=0.2,
            facecolor="auto",
            edgecolor="auto",
            backend=None,
        )

    if args.outFile is None:
        of = sys.stdout
    else:
        of = open(args.outFile, "w")

    df.to_csv(of, sep="\t")

    return 0
