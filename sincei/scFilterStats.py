#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

from deeptools import bamHandler, utilities
from deeptools.mapReduce import mapReduce
from deeptoolsintervals import GTF

import numpy as np
import py2bit
import pandas as pd

# logs
import warnings
import logging

logger = logging.getLogger()
warnings.simplefilter(action="ignore", category=FutureWarning)

from sincei.Utilities import (
    checkMotifs,
    checkGCcontent,
    checkBAMtag,
    checkAlignedFraction,
    getDupFilterTuple,
)
from sincei import ParserCommon


def parseArguments():
    filterParser = ParserCommon.filterOptions()

    io_args = ParserCommon.inputOutputOptions(opts=["bamfiles", "barcodes", "outFile"], requiredOpts=["barcodes"])
    bam_args = ParserCommon.bamOptions(
        suppress_args=["region", "groupTag"],
        default_opts={"binSize": 100000, "distanceBetweenBins": 1000000},
    )
    filter_args = ParserCommon.filterOptions()
    read_args = ParserCommon.readOptions(
        suppress_args=[
            "minFragmentLength",
            "maxFragmentLength",
            "extendReads",
            "centerReads",
        ]
    )
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        parents=[io_args, bam_args, filter_args, read_args, other_args],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
``scFilterStats`` estimates the number of reads that would be filtered given a set of criteria
and prints it to the terminal. Furthermore, it tracks the number of singleton reads.
The following metrics will always be tracked regardless of what you specify (the order output also matches this):

* Total reads (including unmapped)
* Mapped reads
* Reads in blacklisted regions (--blackListFileName)

The following metrics are estimated according to the --binSize and --distanceBetweenBins parameters:

* Estimated mapped reads filtered (the total number of mapped reads filtered for any reason)
* Alignments with a below threshold MAPQ (--minMappingQuality)
* Alignments with at least one missing flag (--samFlagInclude)
* Alignments with undesirable flags (--samFlagExclude)
* Duplicates determined by sincei (--duplicateFilter)
* Duplicates marked externally (e.g., by picard)
* Singletons (paired-end reads with only one mate aligning)
* Wrong strand (due to --filterRNAstrand)

The sum of these may be more than the total number of reads. Note that alignments are sampled from
bins of size --binSize spaced --distanceBetweenBins apart.
""",
        usage="scFilterStats -b sample1.bam sample2.bam -bc barcodes.txt -bl blacklist.bed -o stats.tsv",
        add_help=False,
    )

    return parser


def getFiltered_worker(arglist):
    chrom, start, end, args = arglist
    # Fix the bounds
    if end - start > args.binSize and end - start > args.distanceBetweenBins:
        end -= args.distanceBetweenBins
    if end <= start:
        end = start + 1
    ## open genome if needed
    if args.genome2bit:
        twoBitGenome = py2bit.open(args.genome2bit, True)
    ## open blacklist file
    blackList = None
    if args.blackListFileName is not None:
        blackList = GTF(args.blackListFileName)

    o = []
    for fname in args.bamfiles:
        fh = bamHandler.openBam(fname)
        chromUse = utilities.mungeChromosome(chrom, fh.references)
        prev_pos = set()
        lpos = None

        ## a dict with barcodes = keys
        # metrics
        blacklisted = {}
        minMapq = {}
        samFlagInclude = {}
        samFlagExclude = {}
        internalDupes = {}
        externalDupes = {}
        singletons = {}
        filterRNAstrand = {}
        filterMotifs = {}
        filterGC = {}
        minAlignedFraction = {}
        # trackers
        nFiltered = {}
        total = {}  # This is only used to estimate the percentage affected
        filtered = {}

        for b in args.barcodes:
            total[b] = 0  # This is only used to estimate the percentage affected
            nFiltered[b] = 0
            blacklisted[b] = 0
            minMapq[b] = 0
            samFlagInclude[b] = 0
            samFlagExclude[b] = 0
            internalDupes[b] = 0
            externalDupes[b] = 0
            singletons[b] = 0
            filterRNAstrand[b] = 0
            filterMotifs[b] = 0
            filterGC[b] = 0
            minAlignedFraction[b] = 0

        for read in fh.fetch(chromUse, start, end):
            try:
                bc = read.get_tag(args.cellTag)
            except KeyError:
                continue
            # also keep a counter for barcodes not in whitelist?
            if bc not in args.barcodes:
                continue

            filtered[bc] = 0
            if read.pos < start:
                # ensure that we never double count (in case distanceBetweenBins == 0)
                continue

            if read.flag & 4:
                # Ignore unmapped reads, they were counted already
                continue

            if args.minMappingQuality and read.mapq < args.minMappingQuality:
                filtered[bc] = 1
                minMapq[bc] += 1
            if args.samFlagInclude and read.flag & args.samFlagInclude != args.samFlagInclude:
                filtered[bc] = 1
                samFlagInclude[bc] += 1
            if args.samFlagExclude and read.flag & args.samFlagExclude != 0:
                filtered[bc] = 1
                samFlagExclude[bc] += 1

            if args.minAlignedFraction:
                if not checkAlignedFraction(read, args.minAlignedFraction):
                    filtered[bc] = 1
                    minAlignedFraction[bc] += 1

            ## reads in blacklisted regions
            if blackList and blackList.findOverlaps(
                chrom,
                read.reference_start,
                read.reference_start + read.infer_query_length(always=False) - 1,
            ):
                filtered[bc] = 1
                blacklisted[bc] += 1

            ## Duplicates
            if args.duplicateFilter:
                tup = getDupFilterTuple(read, bc, args.duplicateFilter)
                if lpos is not None and lpos == read.reference_start and tup in prev_pos:
                    filtered[bc] = 1
                    internalDupes[bc] += 1
                if lpos != read.reference_start:
                    prev_pos.clear()
                lpos = read.reference_start
                prev_pos.add(tup)
            if read.is_duplicate:
                filtered[bc] = 1
                externalDupes[bc] += 1
            if read.is_paired and read.mate_is_unmapped:
                filtered[bc] = 1
                singletons[bc] += 1

            ## remove reads with low/high GC content
            if args.GCcontentFilter:
                if not checkGCcontent(read, args.GCcontentFilter[0], args.GCcontentFilter[1]):
                    filtered[bc] = 1
                    filterGC[bc] += 1

            ## remove reads that don't pass the motif filter
            if args.motifFilter:
                test = [checkMotifs(read, chrom, twoBitGenome, m[0], m[1]) for m in args.motifFilter]
                # if none given motif found, return true
                if not any(test):
                    filtered[bc] = 1
                    filterMotifs[bc] += 1

            # filterRNAstrand
            if args.filterRNAstrand:
                if read.is_paired:
                    if args.filterRNAstrand == "forward":
                        if read.flag & 144 == 128 or read.flag & 96 == 64:
                            pass
                        else:
                            filtered[bc] = 1
                            filterRNAstrand[bc] += 1
                    elif args.filterRNAstrand == "reverse":
                        if read.flag & 144 == 144 or read.flag & 96 == 96:
                            pass
                        else:
                            filtered[bc] = 1
                            filterRNAstrand[bc] += 1
                else:
                    if args.filterRNAstrand == "forward":
                        if read.flag & 16 == 16:
                            pass
                        else:
                            filtered[bc] = 1
                            filterRNAstrand[bc] += 1
                    elif args.filterRNAstrand == "reverse":
                        if read.flag & 16 == 0:
                            pass
                        else:
                            filtered[bc] = 1
                            filterRNAstrand[bc] += 1

            total[bc] += 1
            nFiltered[bc] += filtered[bc]
        fh.close()

        # first make a tuple where each entry is a dict of barcodes:value
        tup = (
            total,
            nFiltered,
            blacklisted,
            minMapq,
            samFlagInclude,
            samFlagExclude,
            internalDupes,
            externalDupes,
            singletons,
            filterRNAstrand,
            filterMotifs,
            filterGC,
            minAlignedFraction,
        )

        # now simplify it
        merged = {}
        for b in args.barcodes:
            merged[b] = tuple(merged[b] for merged in tup)
        # now merged is a dict with each key = barcode, values = tuple of stats
        # Now convert it to array
        out = np.stack([v for k, v in merged.items()])
        # out is an array with row = len(barcode) [384], column = len(stats) [11]
        o.append(out)
    return o


def main(args=None):
    args, rowLabels = ParserCommon.validateInputs(parseArguments().parse_args(args))
    if not args.verbose:
        logger.setLevel(logging.CRITICAL)
        warnings.filterwarnings("ignore")

    if args.outFile is None:
        of = sys.stdout
    else:
        of = open(args.outFile, "w")

    for bam in args.bamfiles:
        x = bamHandler.openBam(bam, returnStats=True, nThreads=args.numberOfProcessors)[0]
        chrom_sizes = list(zip(x.references, x.lengths))

        checkBAMtag(x, bam, args.cellTag)
        if args.groupTag:
            checkBAMtag(x, bam, args.groupTag)
            sys.stderr.write("--groupTag is not implemented for scFilterStats yet! \
            Please split your BAM file by {} and re-run scFilterStats. \n".format(args.groupTag))
            exit(1)
        x.close()

    # Get the remaining metrics
    res = mapReduce(
        [args],
        getFiltered_worker,
        chrom_sizes,
        genomeChunkLength=args.binSize + args.distanceBetweenBins,
        blackListFileName=args.blackListFileName,
        numberOfProcessors=args.numberOfProcessors,
        verbose=args.verbose,
    )
    ## res, should be the list of np.arrays of length (len(barcodes) * 9)

    ## final output is an array where nrows = bamfiles*barcodes, ncol = No. of stats
    final_array = np.asarray(res).sum(axis=0)
    ## get final row/col Names (bamnames_barcode)
    colLabels = [
        "Total_sampled",
        "Filtered",
        "Blacklisted",
        "Low_MAPQ",
        "Missing_Flags",
        "Excluded_Flags",
        "Internal_Duplicates",
        "Marked_Duplicates",
        "Singletons",
        "Wrong_strand",
        "Wrong_motif",
        "Unwanted_GC_content",
        "Low_aligned_fraction",
    ]

    final_df = pd.DataFrame(data=np.concatenate(final_array), index=rowLabels, columns=colLabels)
    final_df.index.name = "Cell_ID"
    ## since stats are approximate, present results as %
    final_df.iloc[:, 1:] = final_df.iloc[:, 1:].div(final_df.Total_sampled, axis=0) * 100

    if args.outFile is not None:
        final_df.to_csv(args.outFile, sep="\t")
    else:
        print(final_df.to_string())

    return 0
