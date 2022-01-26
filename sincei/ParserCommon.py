import argparse
import os

def outputOptions(args=None):
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group('Output options')
    group.add_argument('--outFilePrefix', '-o',
                       help='Output file name prefix.',
                       metavar='FILEPREFIX',
                       type=str,
                       required=True)

    return parser

def otherOptions(args=None):
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group('Other options')

    group.add_argument("--help", "-h", action="help",
                          help="show this help message and exit")

    group.add_argument('--verbose', '-v',
                         help='Set to see processing messages.',
                         action='store_true')

#    group.add_argument('--version', action='version',
#                         version='%(prog)s {}'.format(__version__))

    return parser

def plotOptions(args=None):
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group('Plot options')

    group.add_argument('--plotWidth',
                         default=10,
                         type=float,
                         help='Output plot width (in cm). (Default: %(default)s)')

    group.add_argument('--plotHeight',
                         default=10,
                         type=float,
                         help='Output plot height (in cm). (Default: %(default)s)')

    group.add_argument('--plotFileFormat',
                          metavar='FILETYPE',
                          type=str,
                          choices=['png', 'pdf', 'svg', 'eps'],
                          default='png',
                          help='Image format type. If given, this option '
                          'overrides the image format based on the plotFile '
                          'ending. (Default: %(default)s)')
    return parser

def bcOptions(args=None, barcode=True, required=True):
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group('Barcode options')

    if barcode:
        group.add_argument('--barcodes', '-bc',
                               help="A single-column file containing barcodes (whitelist) to be used for the analysis.",
                               metavar="TXT",
                               required=required)
    else:
        group.add_argument('--whitelist', '-w',
                               help="A single-column file containing the whitelist of barcodes to be used",
                               metavar="TXT",
                               default=None,
                               required=required)

    return parser

## Tools: scBulkCoverage, scCountReads, scFilterStats, scJSD
def bamOptions(args=None, binSize=True, blackList=True):
    """
    Arguments for tools that operate on BAM files
    """

    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group('BAM processing options')

    group.add_argument('--tagName', '-tn',
                          metavar='STR',
                          help='Name of the BAM tag from which to extract barcodes.',
                          type=str,
                          default='BC')

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

    group.add_argument('--numberOfProcessors', '-p',
                          help='Number of processors to use. Type "max/2" to '
                          'use half the maximum number of processors or "max" '
                          'to use all available processors. (Default: %(default)s)',
                          metavar="INT",
                          type=numberOfProcessors,
                          default=1,
                          required=False)

    group.add_argument('--region', '-r',
                          help='Region of the genome to limit the operation '
                          'to - this is useful when testing parameters to '
                          'reduce the computing time. The format is '
                          'chr:start:end, for example --region chr10 or '
                          '--region chr10:456700:891000.',
                          metavar="CHR:START:END",
                          required=False,
                          type=genomicRegion)

    if binSize:
        group.add_argument('--binSize', '-bs',
                              help='Size of the bins, in bases, for the output '
                              'of the bigwig/bedgraph file. (Default: %(default)s)',
                              metavar="INT bp",
                              type=int,
                              default=50)

        group.add_argument('--distanceBetweenBins',
                             metavar='INT',
                             help='To reduce the computation time, not every possible genomic '
                             'bin is sampled. This option allows you to set the distance '
                             'between bins actually sampled from. Larger numbers are sufficient '
                             'for high coverage samples, while smaller values are useful for '
                             'lower coverage samples. Note that if you specify a value that '
                             'results in too few (<1000) reads sampled, the value will be '
                             'decreased. (Default: %(default)s)',
                             default=1000000,
                             type=int)
    if blackList:
        group.add_argument('--blackListFileName', '-bl',
                              help="A BED or GTF file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry. Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered. Please note that you should adjust the effective genome size, if relevant.",
                              metavar="BED file",
                              nargs="+",
                              required=False)


    return parser

## Tools: scBulkCoverage, scCountReads
def filterOptions(args=None):
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group('Read Filtering Options')

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
                                'This option checks for the motif at the 5-end of the read and at the 5-overhang in the genome, '
                                'which is useful in identifying reads properly cut by a restriction-enzyme or MNAse. '
                                'For example, if you want to search for an "A" at the 5\'-end of the read and "TA" at 5\'-overhang, '
                                'use "-m \'A,TA\'". Reads not containing the given motif are discarded. ',
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

## Tools: scBulkCoverage, scCountReads
def readOptions():
    """Common arguments related to BAM files and the interpretation
    of the read coverage
    """
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group('Read Processing Options')

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

## helper functions
def process_args(args=None):
    if not args.labels:
        if args.smartLabels:
            args.labels = smartLabels(args.bamfiles)
        else:
            args.labels = args.bamfiles

    if len(args.bamfiles) != len(args.labels):
        sys.exit("The number of labels does not match the number of BAM files.")

    return args

def genomicRegion(string):
    # remove whitespaces using split,join trick
    region = ''.join(string.split())
    if region == '':
        return None
    # remove undesired characters that may be present and
    # replace - by :
    # N.B., the syntax for translate() differs between python 2 and 3
    try:
        region = region.translate(None, ",;|!{}()").replace("-", ":")
    except:
        region = region.translate({ord(i): None for i in ",;|!{}()"})
    if len(region) == 0:
        raise argparse.ArgumentTypeError(
            "{} is not a valid region".format(string))
    return region

def numberOfProcessors(string):
    import multiprocessing
    availProc = multiprocessing.cpu_count()

    if string == "max/2":  # default case
        # by default half of the available processors are used
        numberOfProcessors = int(availProc * 0.5)
    elif string == "max":
        # use all available processors
        numberOfProcessors = availProc
    else:
        try:
            numberOfProcessors = int(string)
        except ValueError:
            raise argparse.ArgumentTypeError(
                "{} is not a valid number of processors".format(string))

        except Exception as e:
            raise argparse.ArgumentTypeError("the given value {} is not valid. "
                                             "Error message: {}\nThe number of "
                                             "available processors in your "
                                             "computer is {}.".format(string, e, availProc))

        if numberOfProcessors > availProc:
            numberOfProcessors = availProc

    return numberOfProcessors
