#!/usr/bin/env python
import argparse
import sys
import os

from deeptools import parserCommon, bamHandler, utilities
from deeptools.mapReduce import mapReduce
from deeptools.utilities import smartLabels
#from deeptools._version import __version__
from deeptoolsintervals import GTF

import numpy as np
import py2bit
import pandas as pd
## own functions
scriptdir=os.path.abspath(os.path.join(__file__, "../../sincei"))
sys.path.append(scriptdir)
from Utilities import *
import ParserCommon

def parseArguments():
    filterParser = ParserCommon.filterOptions()

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
        This tool identifies barcodes present in a BAM files and produces a list. You can optionally filter these barcodes by matching them to a whitelist or
        based on total counts.
        """,
        usage='Example usage: scFilterBarcodes.py -b sample.bam -w whitelist.txt > barcodes_detected.txt')

    required = parser.add_argument_group('Required arguments')
    required.add_argument('--bamfile', '-b',
                          metavar='FILE',
                          help='BAM file',
                          required=True)

    general = parser.add_argument_group('Optional arguments')

    general.add_argument('--outFile', '-o',
                         type=parserCommon.writableFile,
                         help='The file to write results to. By default, results are printed to the console')

    general.add_argument('--whitelist', '-w',
                           help="A single-column file containing the whitelist of barcodes to be used",
                           metavar="TXT",
                           default=None,
                           required=False)

    general.add_argument('--minHammingDist', '-d',
                           help="Minimum hamming distance to match the barcode in whitelist",
                           metavar="INT",
                           type=int,
                           default=0,
                           required=False)

    general.add_argument('--tagName', '-tn',
                          metavar='STR',
                          help='Name of the BAM tag from which to extract barcodes.',
                          type=str,
                          default='BC')

    general.add_argument('--blackListFileName', '-bl',
                           help='A BED or GTF file containing regions that should be excluded from all analyses. '
                           'Currently this works by rejecting genomic chunks that happen to overlap an entry. '
                           'Consequently, for BAM files, if a read partially overlaps a blacklisted region or '
                           'a fragment spans over it, then the read/fragment might still be considered. Please note '
                           'that you should adjust the effective genome size, if relevant.',
                           metavar="BED file",
                           nargs="+",
                           required=False)

    general.add_argument('--minMappingQuality',
                           metavar='INT',
                           help='If set, only reads that have a mapping '
                           'quality score of at least this are '
                           'considered.',
                           type=int)

    general.add_argument('--binSize', '-bs',
                         metavar='INT',
                         help='Length in bases of the window used to count the barcodes. (Default: %(default)s)',
                         default=1000000,
                         type=int)

    general.add_argument('--numberOfProcessors', '-p',
                         help='Number of processors to use. Type "max/2" to '
                         'use half the maximum number of processors or "max" '
                         'to use all available processors. (Default: %(default)s)',
                         metavar="INT",
                         type=parserCommon.numberOfProcessors,
                         default=1,
                         required=False)

    general.add_argument('--verbose', '-v',
                         help='Set to see processing messages.',
                         action='store_true')


    return parser

## get hamming dist between 2 sequences
def ham_dist(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Undefined")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def getFiltered_worker(arglist):
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
        if blackList and blackList.findOverlaps(chrom, read.reference_start, read.reference_start + read.infer_query_length(always=False) - 1):
            continue

        ## get barcode and append to set
        try:
            bc = read.get_tag(args.tagName)
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


def main(args=None):
    args = parseArguments().parse_args(args)

    # open barcode file and read the content in a list
    with open(args.whitelist, 'r') as f:
        barcodes = f.read().splitlines()
    args.whitelist = barcodes

    if args.outFile is None:
        of = sys.stdout
    else:
        of = open(args.outFile, "w")

    bhs = bamHandler.openBam(args.bamfile, returnStats=True, nThreads=args.numberOfProcessors)[0]
    chrom_sizes = list(zip(bhs.references, bhs.lengths))
    bhs.close()

    # Get the remaining metrics
    res = mapReduce([args],
                    getFiltered_worker,
                    chrom_sizes,
                    genomeChunkLength=args.binSize,
                    blackListFileName=args.blackListFileName,
                    numberOfProcessors=args.numberOfProcessors,
                    verbose=args.verbose)
    ## res, should be a list of sets, collapse it to one set
    final_set = list(set().union(*res))

    for x in final_set:
        of.write(str(x)+'\n')

    return 0

if __name__ == "__main__":
    main()
