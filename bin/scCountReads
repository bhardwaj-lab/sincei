#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import numpy as np
from scipy import sparse, io
import re
import pandas as pd
import anndata as ad
import scanpy as sc
from deeptools import parserCommon
from deeptools.utilities import smartLabels
from deeptools._version import __version__

# own functions
scriptdir=os.path.abspath(os.path.join(__file__, "../../sincei"))

sys.path.append(scriptdir)
import ReadCounter as countR
import ParserCommon
old_settings = np.seterr(all='ignore')


def parseArguments(args=None):
    parser = \
        argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description="""
            ``scCountReads`` computes the read coverages per cell barcode for genomic regions in the provided BAM file(s).
            The analysis can be performed for the entire genome by running the program in 'bins' mode.
            If you want to count the read coverage for specific regions only, use the ``features`` mode instead.
            The standard output of ``scCountReads`` is a ".loom" file with counts, along with rowName (features) and colNames (cell barcodes).

            A detailed sub-commands help is available by typing:

              scCountReads bins -h

              scCountReads features -h

          """,
            epilog='example usages:\n'
                   'scCountReads bins --bamfiles file1.bam file2.bam --barcodes whitelist.txt -o results \n\n'
                   'scCountReads features --BED selection.bed --bamfiles file1.bam file2.bam --barcodes whitelist.txt \n'
                   '-o results'
                   ' \n\n',
            conflict_handler='resolve')

    subparsers = parser.add_subparsers(
        title="commands",
        dest='command',
        description='subcommands',
        help='subcommands',
        metavar='')

    read_args = ParserCommon.readOptions(suppress_args=['filterRNAstrand'])
    filter_args = ParserCommon.filterOptions()
    other_args = ParserCommon.otherOptions()

    # bins mode options
    subparsers.add_parser(
        'bins',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[ParserCommon.inputOutputOptions(opts=['bamfiles', 'barcodes', 'outFilePrefix', 'BED'],
                                                 requiredOpts=['barcodes','outFilePrefix'],
                                                 suppress_args=['BED']),
                 parserCommon.gtf_options(suppress=True),
                 ParserCommon.bamOptions(default_opts={'binSize':10000,
                                        'distanceBetweenBins':0}),
                 read_args, filter_args,
                 get_args(), other_args
                 ],
        help="The reads are counted in bins of equal size. The bin size and distance between bins can be adjusted.",
        add_help=False,
        usage='%(prog)s -bs 100000 --bamfiles file1.bam file2.bam -o results \n')

    # BED file arguments
    subparsers.add_parser(
        'features',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[ParserCommon.inputOutputOptions(opts=['bamfiles', 'barcodes', 'outFilePrefix', 'BED'],
                                                 requiredOpts=['barcodes','outFilePrefix', 'BED']),
                 parserCommon.gtf_options(),
                 ParserCommon.bamOptions(suppress_args=['binSize', 'distanceBetweenBins'],
                                        default_opts={'binSize':10000,
                                                      'distanceBetweenBins':0}),
                 read_args, filter_args,
                 get_args(), other_args
                 ],
        help="The user provides a BED/GTF file containing all regions "
             "that should be counted. A common use would be to count scRNA-seq reads on Genes.",
        usage='%(prog)s --BED selection.bed --bamfiles file1.bam file2.bam --barcodes whitelist.txt -o results\n',
        add_help=False)

    return parser


def get_args():
    parser = argparse.ArgumentParser(add_help=False)
    optional = parser.add_argument_group('Misc arguments')
    optional.add_argument('--genomeChunkSize',
                          type=int,
                          default=None,
                          help='Manually specify the size of the genome provided to each processor. '
                          'The default value of None specifies that this is determined by read '
                          'density of the BAM file.')

    optional.add_argument('--outFileFormat',
                          type=str,
                          default='loom',
                          choices=['loom', 'mtx'],
                          help='Output file format. Default is to write an anndata object of name '
                          '<prefix>.loom, which can either be opened in scanpy, or by downstream tools. '
                          '"mtx" refers to the MatrixMarket sparse-matrix format. The output in this case would be '
                          '<prefix>.counts.mtx, along with <prefix>.rownames.txt and <prefix>.colnames.txt')

    return parser


def process_args(args=None):
    args = parseArguments().parse_args(args)

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
            args.motifFilter = [ x.strip(" ").split(",") for x in args.motifFilter ]

    if args.GCcontentFilter:
        gc = args.GCcontentFilter.strip(" ").split(",")
        args.GCcontentFilter = [float(x) for x in gc]

    ## read the barcode file
    with open(args.barcodes, 'r') as f:
        barcodes = f.read().splitlines()
    f.close()

    ## create row/colNames
    if args.outFileFormat == "mtx":
        mtxFile = args.outFilePrefix + ".counts.mtx"
        rowNamesFile = args.outFilePrefix + ".rownames.txt"
        colNamesFile = args.outFilePrefix + ".colnames.txt"
    else:
        rowNamesFile=None

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
        duplicateFilter=args.duplicateFilter,
        center_read=args.centerReads,
        samFlag_include=args.samFlagInclude,
        samFlag_exclude=args.samFlagExclude,
        minFragmentLength=args.minFragmentLength,
        maxFragmentLength=args.maxFragmentLength,
        zerosToNans=False,
        out_file_for_raw_data=rowNamesFile)

    num_reads_per_bin, regionList = c.run(allArgs=args)

    sys.stderr.write("Number of bins "
                     "found: {}\n".format(num_reads_per_bin.shape[0]))

    if num_reads_per_bin.shape[0] < 1:
        exit("ERROR: too few non zero bins found.\n"
             "If using --region please check that this "
             "region is covered by reads.\n")

    ## create colnames (sampleLabel+barcode)
    newlabels = ["{}::{}".format(re.sub('\.bam', '', a), b) for a in args.labels for b in barcodes ]

    ## write mtx/rownames if asked
    if args.outFileFormat == 'mtx':
        f = open(colNamesFile, "w")
        f.write("\n".join(newlabels))
        f.write("\n")
        f.close()
        ## write the matrix as .mtx
        sp = sparse.csr_matrix(num_reads_per_bin)
        io.mmwrite(mtxFile, sp, field="integer")
    else:
        # write anndata
        adata = ad.AnnData(num_reads_per_bin.T)
        adata.obs = pd.DataFrame({"sample":[x.split('::')[-2] for x in newlabels],
                                "barcodes":[x.split('::')[-1] for x in newlabels]
                                },index=newlabels)

        rows=list(regionList)

        adata.var = pd.DataFrame({"chrom":[x.split('_')[0] for x in rows],
                                "start":[x.split('_')[1] for x in rows],
                                "end":[x.split('_')[2] for x in rows]
                                }, index=rows)

        # export as loom
        adata.write_loom(args.outFilePrefix+".loom")

if __name__ == "__main__":
    main()
