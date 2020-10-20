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

def labelOptions(args=None):
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group('Label options')

    group.add_argument('--labels', '-l',
                          metavar='sample1 sample2',
                          help='User defined labels instead of default labels from '
                               'file names. '
                               'Multiple labels have to be separated by a space, e.g. '
                               '--labels sample1 sample2 sample3',
                          nargs='+')

    group.add_argument('--smartLabels',
                          action='store_true',
                          help='Instead of manually specifying labels for the input '
                          'BAM files, this causes sincei to use the file name '
                          'after removing the path and extension.')

    return parser


def filterOptions(args=None):
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group('Filtering Arguments')

    group.add_argument('--tagName', '-tn',
                          metavar='STR',
                          help='Name of the BAM tag from which to extract barcodes.',
                          type=str,
                          default='BC')

    group.add_argument('--duplicateFilter',
                           help='How to filter for duplicates? Different combinations (using start/end/umi) are possible. '
                           'Read start position and read barcode are always considered. Default (None) would consider all reads. '
                           'Note that in case of paired end data, both reads in the fragment are considered (and kept). So if you wish '
                           'to keep only read1, combine this option with samFlagInclude ',
                           type=str,
                           choices=['start_bc', 'start_bc_umi', 'start_end_bc', 'start_end_bc_umi'],
                           default=None)

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


def read_options():
    """Common arguments related to BAM files and the interpretation
    of the read coverage
    """
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group('Read processing options')

    group.add_argument('--extendReads', '-e',
                       help='This parameter allows the extension of reads to '
                       'fragment size. If set, each read is '
                       'extended, without exception.\n'
                       '*NOTE*: This feature is generally NOT recommended for '
                       'spliced-read data, such as RNA-seq, as it would '
                       'extend reads over skipped regions.\n'
                       '*Single-end*: Requires a user specified value for the '
                       'final fragment length. Reads that already exceed this '
                       'fragment length will not be extended.\n'
                       '*Paired-end*: Reads with mates are always extended to '
                       'match the fragment size defined by the two read mates. '
                       'Unmated reads, mate reads that map too far apart '
                       '(>4x fragment length) or even map to different '
                       'chromosomes are treated like single-end reads. The input '
                       'of a fragment length value is optional. If '
                       'no value is specified, it is estimated from the '
                       'data (mean of the fragment size of all mate reads).\n',
                       type=int,
                       nargs='?',
                       const=True,
                       default=False,
                       metavar="INT bp")

    group.add_argument('--centerReads',
                       help='By adding this option, reads are centered with '
                       'respect to the fragment length. For paired-end data, '
                       'the read is centered at the fragment length defined '
                       'by the two ends of the fragment. For single-end data, the '
                       'given fragment length is used. This option is '
                       'useful to get a sharper signal around enriched '
                       'regions.',
                       action='store_true')

    group.add_argument('--minMappingQuality',
                       metavar='INT',
                       help='If set, only reads that have a mapping '
                       'quality score of at least this are '
                       'considered.',
                       type=int,
                       )

    group.add_argument('--samFlagInclude',
                       help='Include reads based on the SAM flag. For example, '
                       'to get only reads that are the first mate, use a flag of 64. '
                       'This is useful to count properly paired reads only once, '
                       'as otherwise the second mate will be also considered for the '
                       'coverage. (Default: %(default)s)',
                       metavar='INT',
                       default=None,
                       type=int,
                       required=False)

    group.add_argument('--samFlagExclude',
                       help='Exclude reads based on the SAM flag. For example, '
                       'to get only reads that map to the forward strand, use '
                       '--samFlagExclude 16, where 16 is the SAM flag for reads '
                       'that map to the reverse strand. (Default: %(default)s)',
                       metavar='INT',
                       default=None,
                       type=int,
                       required=False)

    group.add_argument('--minFragmentLength',
                       help='The minimum fragment length needed for read/pair '
                       'inclusion. This option is primarily useful '
                       'in ATACseq experiments, for filtering mono- or '
                       'di-nucleosome fragments. (Default: %(default)s)',
                       metavar='INT',
                       default=0,
                       type=int,
                       required=False)

    group.add_argument('--maxFragmentLength',
                       help='The maximum fragment length needed for read/pair '
                       'inclusion. (Default: %(default)s)',
                       metavar='INT',
                       default=0,
                       type=int,
                       required=False)

    return parser


def process_args(args=None):

    args = parse_arguments().parse_args(args)

    if not args.labels:
        if args.smartLabels:
            args.labels = smartLabels(args.bamfiles)
        else:
            args.labels = args.bamfiles

    if len(args.bamfiles) != len(args.labels):
        sys.exit("The number of labels does not match the number of BAM files.")

    return args
