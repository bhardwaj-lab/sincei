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
sys.path.append('/home/vbhardwaj/programs/sincei/sincei')
from utilities import *

def getStats_worker(arglist):
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
        ## initiate a dict with all values to keep per read
        info_list = []#dict.fromkeys(['barcode', 'position', 'duplicate', 'GCcontent', 'strand'])

        for read in fh.fetch(chromUse, start, end):
            ## general filtering
            if read.pos < start:
                # ensure that we never double count (in case distanceBetweenBins == 0)
                continue
            if read.flag & 4:
                # Ignore unmapped reads, they were counted already
                continue
            if args.minMappingQuality and read.mapq < args.minMappingQuality:
                continue
            if args.minAlignedFraction:
                if not checkAlignedFraction(read, args.minAlignedFraction):
                    continue
            if blackList and blackList.findOverlaps(chrom, read.reference_start, read.reference_start + read.infer_query_length(always=False) - 1):
                continue
            if args.motifFilter:
                test = [ checkMotifs(read, chrom, twoBitGenome, m[0], m[1]) for m in args.motifFilter ]
                # if none given motif found, return true
                if not any(test):
                    continue

            # now collect info
            info = [None for n in range(0,8)]
            info[0] = chromUse
            if args.barcodes is not None:
                bc = read.get_tag(args.tagName)
                if bc in args.barcodes:
                    info[1] = True
                else:
                    info[1] = False
            else:
                info[1] = False
                bc = None
            info[2] = read.reference_start

            ## Duplicates
            tup = getDupFilterTuple(read, bc, args.duplicateFilter)
            if lpos is not None and lpos == read.reference_start \
                    and tup in prev_pos:
                info[3] = True# read is duplicate
                info[4] = 0
            else:
                info[3] = False
            if lpos != read.reference_start:
                prev_pos.clear()
                ## add distance to last read
                if lpos is None:
                    info[4] = 0
                else:
                    info[4] = abs(float(lpos - read.reference_start))
                ## also add the distance of this read to previous pos

            lpos = read.reference_start
            prev_pos.add(tup)

            info[5] = checkGCcontent(read, args.GCcontentFilter[0], args.GCcontentFilter[1], returnGC=True)
            # filterRNAstrand
            info[6] = read.is_reverse
            if args.getReadID:
                info[7] = read.query_name
            info_list.append(info)
        fh.close()

        # out is an array with row = len(barcode) [384], column = len(stats) [11]
        o.append(info_list)
    return o