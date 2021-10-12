#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import numpy as np
from scipy import sparse, io

import pandas as pd
import anndata as ad

from deeptools import parserCommon
from deeptools.utilities import smartLabels
from deeptools._version import __version__

# own functions
scriptdir=os.path.join(os.path.dirname(os.path.dirname(__file__)), "sincei")
print(scriptdir)
sys.path.append(scriptdir)
import ReadCounter as countR
import ParserCommon

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
    read_options_parser = ParserCommon.read_options()

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
    outputParser = ParserCommon.output()
    filterParser = ParserCommon.filterOptions()
    label_parser = ParserCommon.labelOptions()
    parser = argparse.ArgumentParser(parents=[outputParser, filterParser, label_parser],
                                    add_help=False)
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

    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument("--help", "-h", action="help",
                          help="show this help message and exit")

    optional.add_argument('--genomeChunkSize',
                          type=int,
                          default=None,
                          help='Manually specify the size of the genome provided to each processor. '
                          'The default value of None specifies that this is determined by read '
                          'density of the BAM file.')

    optional.add_argument('--outFileFormat',
                          type=str,
                          default='anndata',
                          choices=['anndata', 'mtx'],
                          help='Output file format. Default is to write an anndata object of name '
                          '<prefix>.h5ad, which can either be opened in scanpy, or by downstream tools. '
                          '"mtx" refers to the MatrixMarket sparse-matrix format. The output in this case would be '
                          '<prefix>.counts.mtx, along with <prefix>.rownames.txt and <prefix>.colnames.txt')

    if case == 'bins':
        optional.add_argument('--binSize', '-bs',
                              metavar='INT',
                              help='Length in bases of the window used '
                                   'to sample the genome. (Default: %(default)s)',
                              default=10000,
                              type=int)

        optional.add_argument('--distanceBetweenBins', '-n',
                              metavar='INT',
                              help='By default, scCountReads considers consecutive '
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
    newlabels = ["{}_{}".format(a, b) for a in args.labels for b in barcodes ]

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
        adata.obs = pd.DataFrame({"sample":[x.split('_')[0] for x in newlabels],
                                "barcodes":[x.split('_')[1] for x in newlabels]
                                },index=newlabels)

        rows=list(regionList)

        adata.var = pd.DataFrame({"chrom":[x.split('_')[0] for x in rows],
                                "start":[x.split('_')[1] for x in rows],
                                "end":[x.split('_')[2] for x in rows]
                                }, index=rows)

        ## add QC stats to the anndata object
        # 1. scanpy metrics
        sc.pp.calculate_qc_metrics(adata, inplace=True)

        # 2. fraction of regions and cells with signal
        # how many regions have non-zero counts per cell, as % of total regions
        nonzero_cells_per_region=np.sum(adata.X > 0, axis=0)
        nonzero_frac_region = [float(x)/adata.obs.shape[0] for x in nonzero_cells_per_region.tolist()[0]]
        adata.var['fraction_cells_with_signal'] = nonzero_frac_region
        # how many regions have non-zero counts per cell, as % of total regions
        nonzero_regions_per_cell=np.sum(adata.X > 0, axis=1)
        nonzero_frac = [float(x)/adata.var.shape[0] for x in nonzero_regions_per_cell]
        chic_adata.obs['fraction_regions_with_signal'] = nonzero_frac

        # 3. Gini coefficient
        gini_list=[]
        for i in range(adata.shape[0]):
            ar=adata.X[:,i].todense()
            ar=np.array(ar).flatten()
            if len(ar[ar>0]) > 2:
                gini_list.append(gini(ar[ar>0]))
            else:
                gini_list.append(1.0)
        adata.obs['gini_coefficient'] = li

        # export as loom
        adata.write_h5ad(args.outFilePrefix+".h5ad")

if __name__ == "__main__":
    main()
