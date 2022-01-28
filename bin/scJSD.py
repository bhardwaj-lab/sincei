#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import time
import sys
import argparse
import numpy as np
import pandas as pd
from deeptools.plotFingerprint import getSyntheticJSD
from deeptools import parserCommon

from matplotlib.pyplot import plot
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

scriptdir=os.path.abspath(os.path.join(__file__, "../../sincei"))
sys.path.append(scriptdir)
## own functions
import ReadCounter as countR
import ParserCommon

## plot KDE of JSD values
from numpy import array, linspace
from sklearn import neighbors

old_settings = np.seterr(all='ignore')
MAXLEN = 10000000


def get_args():
    parser = argparse.ArgumentParser(add_help=False,
                                     conflict_handler='resolve')
    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument('--numberOfSamples', '-n',
                          help='The number of bins that are sampled from the genome, '
                          'for which the overlapping number of reads is computed. (Default: %(default)s)',
                          default=1e5,
                          type=int)

    optional.add_argument('--skipZeros',
                          help='If set, then regions with zero overlapping reads'
                          'for *all* given BAM files are ignored. This '
                          'will result in a reduced number of read '
                          'counts than that specified in --numberOfSamples',
                          action='store_true')


    return parser

def parse_arguments(args=None):
    io_args = ParserCommon.inputOutputOptions(opts=['bamfiles', 'barcodes', 'outFile'],
                                             requiredOpts=['bamfiles', 'barcodes', 'outFile'])
    bam_args = ParserCommon.bamOptions(suppress_args=['region', 'distanceBetweenBins'],
                                      default_opts={'binSize': 10000})
    read_args = ParserCommon.readOptions(suppress_args=['filterRNAstrand', 'extendReads', 'centerReads'])
    filter_args = ParserCommon.filterOptions(suppress_args=['motifFilter', 'genome2bit', 'GCcontentFilter', 'minAlignedFraction'])
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        parents=[io_args, bam_args, read_args,
                 filter_args, get_args(), other_args],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='This tool samples regions in the genome from BAM files '
        'and compares the cumulative read coverages for each cell on those regions. '
        'to a synthetic cell with poisson distributed reads using Jansson Shannon Distance. '
        'Cells with high enrichment of signals show a higher JSD compared to cells whose signal '
        'is homogenously distrubuted.',
        conflict_handler='resolve',
        usage='An example usage is: plotFingerprint -b treatment.bam control.bam '
        '-plot fingerprint.png',
        add_help=False)

    return parser


def main(args=None):
    args = ParserCommon.process_args(parse_arguments().parse_args(args))

    ## read the barcode file
    with open(args.barcodes, 'r') as f:
        barcodes = f.read().splitlines()
    f.close()

    ## Count
    c = countR.CountReadsPerBin(
        args.bamfiles,
        binLength=args.binSize,
        numberOfSamples=args.numberOfSamples,
        barcodes=barcodes,
        tagName=args.tagName,
        motifFilter=None,
        genome2bit=None,
        GCcontentFilter=None,
        blackListFileName=args.blackListFileName,
        numberOfProcessors=args.numberOfProcessors,
        verbose=args.verbose,
        region=None,
        bedFile=None,
        extendReads=args.extendReads,
        minMappingQuality=args.minMappingQuality,
        duplicateFilter=args.duplicateFilter,
        center_read=False,
        samFlag_include=args.samFlagInclude,
        samFlag_exclude=args.samFlagExclude,
        minFragmentLength=args.minFragmentLength,
        maxFragmentLength=args.maxFragmentLength,
        zerosToNans=False,
        sumCoveragePerBin=True)

    num_reads_per_bin, _ = c.run(allArgs=None)

    if num_reads_per_bin.sum() == 0:
        import sys
        sys.stderr.write(
            "\nNo reads were found in {} regions sampled. Check that the\n"
            "min mapping quality is not overly high and that the \n"
            "chromosome names between bam files are consistant.\n"
            "For small genomes, decrease the --numberOfSamples.\n"
            "\n".format(num_reads_per_bin.shape[0]))
        exit(1)

    if args.skipZeros:
        num_reads_per_bin = countR.remove_row_of_zeros(num_reads_per_bin)

    total = len(num_reads_per_bin[:, 0])
    x = np.arange(total).astype('float') / total  # normalize from 0 to 1

    jsd_all = []
    for i in range(0, num_reads_per_bin.shape[1]):
        jsd_all.append(getSyntheticJSD(num_reads_per_bin[:, i]))

    ## create colnames (sampleLabel+barcode)
    newlabels = ["{}_{}".format(a, b) for a in args.labels for b in barcodes ]
    df = pd.DataFrame({'cell': newlabels, 'jsd': jsd_all})
    df.to_csv(args.outFile, sep = "\t")



if __name__ == "__main__":
    main()
