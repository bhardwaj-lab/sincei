import argparse
import os
#from scDeepTools._version import __version__

def output(args=None):
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group('Output options')
    group.add_argument('--outFilePrefix', '-o',
                       help='Output file name prefix.',
                       metavar='FILEPREFIX',
                       type=str,
                       required=True)

    return parser

def filterOptions(args=None):
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group('Filtering Arguments')

    group.add_argument('--tagName', '-tn',
                          metavar='STR',
                          help='Name of the BAM tag from which to extract barcodes.',
                          type=str,
                          default='BC')

    group.add_argument('--motifFilter', '-m',
                          metavar='STR',
                          help='Check whether a given motif is present in the read and the corresponding reference genome. '
                                'This function checks for the motif at the 5-end of the read and at the 5-overhang in the genome. '
                                'Reads not containing the given motif are not discarded. ',
                          type=str,
                          nargs='+',
                          default=None)

    group.add_argument('--genome2bit', '-g',
                          metavar='STR',
                          help='If --motifFilter is provided, please also provide the genome sequence (in 2bit format). ',
                          type=str,
                          default=None)

    group.add_argument('--GCcontentFilter', '-gc',
                          metavar='STR',
                          help='Check whether the GC content of the read falls within the provided range. '
                                'If the GC content of the reads fall outside the range, they are discarded. ',
                          type=str,
                          default=None)

    group.add_argument('--minAlignedFraction',
                           help='Minimum fraction of the reads which should be aligned to be counted. This includes '
                           'mismatches tolerated by the aligners, but excludes InDels/Clippings (Default: %(default)s)',
                           metavar='FLOAT',
                           default=None,
                           type=float,
                           required=False)
    return parser
