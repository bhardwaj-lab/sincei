#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import numpy as np
from scipy import sparse, io

from deeptools import parserCommon
from deeptools.utilities import smartLabels
from deeptools._version import __version__

# own functions
import scReadCounter as countR

old_settings = np.seterr(all='ignore')


def parse_arguments(args=None):
    parser = \
        argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description="""

``scCountReads`` computes the read coverages per cell barcode for genomic regions in the provided BAM file(s).
The analysis can be performed for the entire genome by running the program in 'bins' mode.
If you want to count the read coverage for specific regions only, use the ``BED-file`` mode instead.
The standard output of ``scCountReads`` is a ".mtx" file with counts, along with rowName and colNames in a single-column .txt file.

A detailed sub-commands help is available by typing:

  scCountReads bins -h

  scCountReads BED-file -h

""",
            epilog='example usages:\n'
                   'scCountReads bins --bamfiles file1.bam file2.bam --barcodes whitelist.txt -o results \n\n'
                   'scCountReads BED-file --BED selection.bed --bamfiles file1.bam file2.bam --barcodes whitelist.txt \n'
                   '-o results'
                   ' \n\n',
            conflict_handler='resolve')

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))
    subparsers = parser.add_subparsers(
        title="commands",
        dest='command',
        description='subcommands',
        help='subcommands',
        metavar='')

    parent_parser = parserCommon.getParentArgParse(binSize=False)
    read_options_parser = parserCommon.read_options()

    # bins mode options
    subparsers.add_parser(
        'bins',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[bamcorrelate_args(case='bins'),
                 parent_parser, read_options_parser,
                 parserCommon.gtf_options(suppress=True)
                 ],
        help="The coverage calculation is done for consecutive bins of equal "
             "size (10 kilobases by default). The bin size and distance between bins can be adjusted.",
        add_help=False,
        usage='%(prog)s '
              '--bamfiles file1.bam file2.bam '
              '-o results.npz \n')

    # BED file arguments
    subparsers.add_parser(
        'BED-file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[bamcorrelate_args(case='BED-file'),
                 parent_parser, read_options_parser,
                 parserCommon.gtf_options()
                 ],
        help="The user provides a BED/GTF file containing all regions "
             "that should be considered for the coverage analysis. A "
             "common use would be to count scRNA-seq coverage on Genes.",
        usage='%(prog)s --BED selection.bed --bamfiles file1.bam file2.bam --barcodes whitelist.txt -o results\n',
        add_help=False)

    return parser


def bamcorrelate_args(case='bins'):
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group('Required arguments')

    # define the arguments
    required.add_argument('--bamfiles', '-b',
                          metavar='FILE1 FILE2',
                          help='List of indexed bam files separated by spaces.',
                          nargs='+',
                          required=True)

    required.add_argument('--barcodes', '-bc',
                          metavar='TXT',
                          help='A single-column text file with barcodes (whitelist) to count.',
                          required=True)

    required.add_argument('--outFilePrefix', '-out', '-o',
                          help='Prefix of file name to save the results. The output is a sparse matrix '
                               'which can be subsequently loaded into R or python for further analysis.',
                          type=str)

    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument("--help", "-h", action="help",
                          help="show this help message and exit")

    optional.add_argument('--tagName', '-tn',
                          metavar='STR',
                          help='Name of the BAM tag from which to extract barcodes.',
                          type=str,
                          default='BC')

    optional.add_argument('--motifFilter', '-m',
                          metavar='STR',
                          help='Check whether a given motif is present in the read and the corresponding reference genome. '
                                'This function checks for the motif at the 5-end of the read and at the 5-overhang in the genome. '
                                'Reads not containing the given motif are not discarded. ',
                          type=str,
                          default=None)

    optional.add_argument('--genome2bit', '-g',
                          metavar='STR',
                          help='If --motifFilter is provided, please also provide the genome sequence (in 2bit format). ',
                          type=str,
                          default=None)

    optional.add_argument('--GCcontentFilter', '-gc',
                          metavar='STR',
                          help='Check whether the GC content of the read falls within the provided range. '
                                'If the GC content of the reads fall outside the range, they are discarded. ',
                          type=str,
                          default=None)

    optional.add_argument('--labels', '-l',
                          metavar='sample1 sample2',
                          help='User defined labels instead of default labels from '
                               'file names. '
                               'Multiple labels have to be separated by a space, e.g. '
                               '--labels sample1 sample2 sample3',
                          nargs='+')

    optional.add_argument('--smartLabels',
                          action='store_true',
                          help='Instead of manually specifying labels for the input '
                          'BAM files, this causes deepTools to use the file name '
                          'after removing the path and extension.')

    optional.add_argument('--genomeChunkSize',
                          type=int,
                          default=None,
                          help='Manually specify the size of the genome provided to each processor. '
                          'The default value of None specifies that this is determined by read '
                          'density of the BAM file.')

    if case == 'bins':
        optional.add_argument('--binSize', '-bs',
                              metavar='INT',
                              help='Length in bases of the window used '
                                   'to sample the genome. (Default: %(default)s)',
                              default=10000,
                              type=int)

        optional.add_argument('--distanceBetweenBins', '-n',
                              metavar='INT',
                              help='By default, multiBamSummary considers consecutive '
                              'bins of the specified --binSize. However, to '
                              'reduce the computation time, a larger distance '
                              'between bins can by given. Larger distances '
                              'result in fewer bins considered. (Default: %(default)s)',
                              default=0,
                              type=int)

        required.add_argument('--BED',
                              help=argparse.SUPPRESS,
                              default=None)
    else:
        optional.add_argument('--binSize', '-bs',
                              help=argparse.SUPPRESS,
                              default=10000,
                              type=int)

        optional.add_argument('--distanceBetweenBins', '-n',
                              help=argparse.SUPPRESS,
                              metavar='INT',
                              default=0,
                              type=int)

        required.add_argument('--BED',
                              help='Limits the coverage analysis to '
                              'the regions specified in these files.',
                              metavar='FILE1.bed FILE2.bed',
                              nargs='+',
                              required=True)

    group = parser.add_argument_group('Output optional options')

    return parser


def process_args(args=None):
    args = parse_arguments().parse_args(args)

    if args.labels and len(args.bamfiles) != len(args.labels):
        print("The number of labels does not match the number of bam files.")
        exit(0)
    if not args.labels:
        if args.smartLabels:
            args.labels = smartLabels(args.bamfiles)
        else:
            args.labels = [os.path.basename(x) for x in args.bamfiles]

    return args


def main(args=None):
    """
    1. get read counts at different positions either
    all of same length or from genomic regions from the BED file

    2. save data for further plotting

    """
    args = process_args(args)

    if 'BED' in args:
        bed_regions = args.BED
    else:
        bed_regions = None

    ## Motif and GC filter
    if args.motifFilter:
        if not args.genome2bit:
            print("MotifFilter asked but genome (2bit) file not provided.")
            sys.exit(1)
        else:
            args.motifFilter = args.motifFilter.strip(" ").split(",")

    if args.GCcontentFilter:
        gc = args.GCcontentFilter.strip(" ").split(",")
        args.GCcontentFilter = [float(x) for x in gc]

    ## read the barcode file
    with open(args.barcodes, 'r') as f:
        barcodes = f.read().splitlines()
    f.close()

    ## create row/colNames
    mtxFile = args.outFilePrefix + ".counts.mtx"
    rowNamesFile = args.outFilePrefix + ".rownames.txt"
    colNamesFile = args.outFilePrefix + ".colnames.txt"

    stepSize = args.binSize + args.distanceBetweenBins
    c = countR.CountReadsPerBin(
        args.bamfiles,
        binLength=args.binSize,
        stepSize=stepSize,
        barcodes=barcodes,
        tagName=args.tagName,
        motifFilter=args.motifFilter,
        genome2bit=args.genome2bit,
        GCcontentFilter=args.GCcontentFilter,
        numberOfSamples=None,
        genomeChunkSize=args.genomeChunkSize,
        numberOfProcessors=args.numberOfProcessors,
        verbose=args.verbose,
        region=args.region,
        bedFile=bed_regions,
        blackListFileName=args.blackListFileName,
        extendReads=args.extendReads,
        minMappingQuality=args.minMappingQuality,
        ignoreDuplicates=args.ignoreDuplicates,
        center_read=args.centerReads,
        samFlag_include=args.samFlagInclude,
        samFlag_exclude=args.samFlagExclude,
        minFragmentLength=args.minFragmentLength,
        maxFragmentLength=args.maxFragmentLength,
        zerosToNans=False,
        out_file_for_raw_data=rowNamesFile)

    num_reads_per_bin = c.run(allArgs=args)

    sys.stderr.write("Number of bins "
                     "found: {}\n".format(num_reads_per_bin.shape[0]))

    #if num_reads_per_bin.shape[0] < 2:
    #    exit("ERROR: too few non zero bins found.\n"
    #         "If using --region please check that this "
    #         "region is covered by reads.\n")

    ## add labels to barcodes and write
    newlabels = ["{}_{}".format(a, b) for a in args.labels for b in barcodes ]

    f = open(colNamesFile, "w")
    f.write("\n".join(newlabels))
    f.write("\n")
    f.close()

    ## write the matrix as .mtx
    sp = sparse.csr_matrix(num_reads_per_bin)
    io.mmwrite(mtxFile, sp, field="integer")


if __name__ == "__main__":
    main()
