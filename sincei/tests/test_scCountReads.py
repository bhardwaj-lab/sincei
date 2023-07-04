from sincei.scCountReads import *
from sincei import ReadCounter as countR
from sincei.Utilities import *

import pandas as pd
import numpy as np
import numpy.testing as nt

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/_data/"

def getCountReadsArgs(T):
    "return arg object based on input T"
    if T == "bins":
        TYPE="bins -bs 10000"
    elif T == "gtf":
        TYPE="features --BED {0}/Ogfrl.gtf".format(ROOT)
    elif T == "bed":
        TYPE="features --BED {0}/test_regions.bed".format(ROOT)

    Args = "{1} -b {0}/SL2-1.bam {0}/SL2-2.bam" \
            " -bc {0}/test_barcodes.txt -ct BC -o None --smartLabels " \
            "--region chr1:23365000:23385000".format(ROOT, TYPE).split()
    args, newlabels = ParserCommon.validateInputs(parseArguments().parse_args(Args))

    return args, newlabels

def getCountReadsOutput(arg, dedup):
    """
    Setup the CR object based on testdata (2x bams, 5 barcodes) and input args
    """
    args, newlabels = getCountReadsArgs(arg)

    stepSize = args.binSize + args.distanceBetweenBins
    args.duplicateFilter = dedup

    c = countR.CountReadsPerBin(
        args.bamfiles,
        binLength=args.binSize,
        stepSize=stepSize,
        barcodes=args.barcodes,
        cellTag=args.cellTag,
        groupTag=args.groupTag,
        groupLabels=newlabels,
        numberOfSamples=None,
        region=args.region,
        bedFile=args.BED,
        duplicateFilter=args.duplicateFilter,
        zerosToNans=False,
        out_file_for_raw_data=None,
    )
    num_reads_per_bin, regionList = c.run(allArgs=args)

    return num_reads_per_bin, regionList


def getExpectedOutput(T, dedup):
    if T == 'bins':
        regions = np.array(['chr1_23360000_23370000::None',
                                  'chr1_23370000_23380000::None',
                                  'chr1_23380000_23385000::None'])
        if dedup is None:
            counts = np.array([
                [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                [14.,  0.,  0., 32., 10.,  2.,  0., 22.,  4.,  3.],
                [ 0.,  6.,  4.,  0.,  0.,  0.,  6.,  4.,  0.,  8.]
                ])
        if dedup == 'start_bc_umi':
            counts = np.array([
                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                [6., 0., 0., 8., 6., 2., 0., 5., 2., 2.],
                [0., 3., 2., 0., 0., 0., 3., 4., 0., 5.]
                ])

    elif T == 'bed':
        regions = np.array(['chr1_23365000_23377000::chr1:23365000-23377000',
                                  'chr1_23365000_23377000::chr1:23365000-23377000_r1',
                                  'chr1_23365000_23385000::chr1:23365000-23385000'])
        if dedup is None:
            counts = np.array([
                [ 0.,  0.,  0., 32., 10.,  2.,  0.,  0.,  0.,  0.],
                [ 0.,  0.,  0., 32., 10.,  2.,  0.,  0.,  0.,  0.],
                [14.,  6.,  4., 32., 10.,  2.,  6., 26.,  4., 10.]
                ])
        if dedup == 'start_bc_umi':
            counts = np.array([
                [0., 0., 0., 8., 6., 2., 0., 0., 0., 0.],
                [0., 0., 0., 8., 6., 2., 0., 0., 0., 0.],
                [6., 3., 2., 8., 6., 2., 3., 9., 2., 6.]
            ])

    elif T == 'gtf':
        regions = np.array(['chr1_23366423_23383175::ENSMUST00000027343.5',
                                  'chr1_23369723_23380801::ENSMUST00000186064.6',
                                  'chr1_23370267_23397541::ENSMUST00000188677.1'])
        if dedup is None:
            counts = np.array([
                [14.,  6.,  0., 32., 10.,  2.,  6., 26.,  4., 10.],
                [14.,  6.,  0., 32., 10.,  2.,  0., 26.,  4., 10.],
                [14.,  6.,  4., 32.,  4.,  2.,  6., 26.,  4., 10.]
                ])
        if dedup == 'start_bc_umi':
             counts = np.array([
                [6., 3., 0., 8., 6., 2., 3., 9., 2., 6.],
                [6., 3., 0., 8., 6., 2., 0., 9., 2., 6.],
                [6., 3., 2., 8., 4., 2., 3., 9., 2., 6.]
             ])

    return counts, regions


def testCountReads_bins_default():
    # Expected output
    valid_counts, valid_regions = getExpectedOutput('bins', None)
    # Actual output
    observed_counts, observed_regions = getCountReadsOutput('bins', None)
    # Test
    nt.assert_array_equal(valid_regions, observed_regions)
    nt.assert_array_equal(valid_counts, observed_counts)

def testCountReads_bins_dedup():
    # Expected output
    valid_counts, valid_regions = getExpectedOutput('bins', 'start_bc_umi')
    # Actual output
    observed_counts, observed_regions = getCountReadsOutput('bins', 'start_bc_umi')
    # Test
    nt.assert_array_equal(valid_regions, observed_regions)
    nt.assert_array_equal(valid_counts, observed_counts)

def testCountReads_bed_default():
    # Expected output
    valid_counts, valid_regions = getExpectedOutput('bed', None)
    # Actual output
    observed_counts, observed_regions = getCountReadsOutput('bed', None)
    # Test
    nt.assert_array_equal(valid_regions, observed_regions)
    nt.assert_array_equal(valid_counts, observed_counts)

def testCountReads_bed_dedup():
    # Expected output
    valid_counts, valid_regions = getExpectedOutput('bed', 'start_bc_umi')
    # Actual output
    observed_counts, observed_regions = getCountReadsOutput('bed', 'start_bc_umi')
    # Test
    nt.assert_array_equal(valid_regions, observed_regions)
    nt.assert_array_equal(valid_counts, observed_counts)

def testCountReads_gtf_default():
    # Expected output
    valid_counts, valid_regions = getExpectedOutput('gtf', None)
    # Actual output
    observed_counts, observed_regions = getCountReadsOutput('gtf', None)
    # Test
    nt.assert_array_equal(valid_regions, observed_regions)
    nt.assert_array_equal(valid_counts, observed_counts)

def testCountReads_gtf_dedup():
    # Expected output
    valid_counts, valid_regions = getExpectedOutput('gtf', 'start_bc_umi')
    # Actual output
    observed_counts, observed_regions = getCountReadsOutput('gtf', 'start_bc_umi')
    # Test
    nt.assert_array_equal(valid_regions, observed_regions)
    nt.assert_array_equal(valid_counts, observed_counts)
