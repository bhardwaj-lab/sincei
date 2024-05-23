#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import os
import itertools
import numpy as np
import pandas as pd
from deeptools.getScaleFactor import get_scale_factor
from deeptools.bamHandler import openBam

# logs
import warnings
import logging

logger = logging.getLogger()
warnings.simplefilter(action="ignore", category=FutureWarning)

## own Functions
# scriptdir=os.path.abspath(os.path.join(__file__, "../../sincei"))
# sys.path.append(scriptdir)

from sincei import ParserCommon
from sincei import WriteBedGraph

debug = 0


def parseArguments():
    io_args = ParserCommon.inputOutputOptions(opts=["bamfiles", "groupInfo", "outFilePrefix"])
    bam_args = ParserCommon.bamOptions(default_opts={"binSize": 100})
    read_args = ParserCommon.readOptions()
    filter_args = ParserCommon.filterOptions()
    other_args = ParserCommon.otherOptions()
    parser = argparse.ArgumentParser(
        parents=[io_args, get_args(), bam_args, filter_args, read_args, other_args],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This tool takes alignments of reads or fragments "
        "as input (BAM files), along with cell grouping information, such as "
        "barcode -> batch, or barcode -> cluster, as tsv file, and generates a "
        "coverage track (bigWig or bedGraph) per group as output. "
        "The coverage is calculated as the number of reads per bin, "
        "where bins are short consecutive counting windows of a defined "
        "size. It is possible to extended/change the length of the reads "
        "to better reflect the actual fragment length. *scBulkCoverage* "
        "offers normalization per cluster using different methods \n",
        usage="An example usage is:"
        "$ scBulkCoverage -b file1.bam file2.bam --labels file1 file2 -g scClusterCells_output.tsv -o coverage.bw",
        add_help=False,
    )

    return parser


def get_args():
    parser = argparse.ArgumentParser(add_help=False)

    optional = parser.add_argument_group("Coverage-related Options")
    optional.add_argument(
        "--outFileFormat",
        "-of",
        help='Output file type. Either "bigwig" or "bedgraph".',
        choices=["bigwig", "bedgraph"],
        default="bigwig",
    )

    optional.add_argument(
        "--normalizeUsing",
        "-n",
        help="How to normalize the pseudo-bulk counts. Options are "
        '"CPM": normalized each bin to the counts per million mapped reads in that group.\n'
        '"Frequency": binarize the coverage per bin and normalize to the total no. of cells per group. \n'
        '"Mean": get mean signal per bin across cells in each group.\n'
        '"None": simply return the sum of coverage per group.',
        choices=["CPM", "Frequency", "Mean", "None"],
        default="CPM",
    )

    optional.add_argument(
        "--ignoreForNormalization",
        "-ig",
        help="Chromosomes to skip while calculating normalization factors",
        nargs="+",
        default=None,
    )

    optional.add_argument(
        "--normalizeByReference",
        "-nr",
        help="NOT IMPLEMENTED: Normalize each group of cells by a reference group (which must be present in the --groupinfo file)"
        "Note that the --normalizeUsing method is applied beforehand.",
        choices=["ratio", "log2_ratio", "difference", "None"],
        default=None,
    )

    optional.add_argument(
        "--scaleFactor",
        help="The computed scaling factor (or 1, if not applicable) will "
        "be multiplied by this. (Default: %(default)s)",
        default=1.0,
        type=float,
        required=False,
    )

    optional.add_argument(
        "--MNase",
        help="Determine nucleosome positions from MNase-seq/CUTnRUN data. "
        "Only 3 nucleotides at the center of each fragment are counted. "
        "The fragment ends are defined by the two mate reads. Only fragment lengths"
        "between 130 - 200 bp are considered to avoid dinucleosomes or other artifacts. "
        "By default, any fragments smaller or larger than this are ignored. To "
        "over-ride this, use the --minFragmentLength and --maxFragmentLength options, "
        "which will default to 130 and 200 if not otherwise specified in the presence "
        "of --MNase. *NOTE*: Requires paired-end data. A bin size of 1 is recommended.",
        action="store_true",
    )

    optional.add_argument(
        "--Offset",
        help="Uses this offset inside of each read as the signal. This is useful in "
        "cases like RiboSeq or GROseq, where the signal is 12, 15 or 0 bases past the "
        "start of the read. This can be paired with the --filterRNAstrand option. "
        "Note that negative values indicate offsets from the end of each read. A value "
        "of 1 indicates the first base of the alignment (taking alignment orientation "
        "into account). Likewise, a value of -1 is the last base of the alignment. An "
        "offset of 0 is not permitted. If two values are specified, then they will be "
        "used to specify a range of positions. Note that specifying something like "
        "--Offset 5 -1 will result in the 5th through last position being used, which "
        "is equivalent to trimming 4 bases from the 5-prime end of alignments. Note "
        "that if you specify --centerReads, the centering will be performed before the "
        "offset.",
        metavar="INT",
        type=int,
        nargs="+",
        required=False,
    )

    return parser


def main(args=None):
    args = ParserCommon.validateInputs(
        parseArguments().parse_args(args), process_barcodes=False
    )  # newlabels are discarded, since they are re-created from groupInfo file
    global debug
    if args.verbose:
        sys.stderr.write("Specified --scaleFactor: {}\n".format(args.scaleFactor))
        debug = 1
    else:
        debug = 0
        logger.setLevel(logging.CRITICAL)
        warnings.filterwarnings("ignore")

    ## read the group info file (TODO: write a validator)
    df_err = """ *Error*:No. of columns in --groupInfo file not recognized.
             Please provide either 3 (sample, barcode, group) or 4 (sample::barcode, umap1, umap2, group) column file"""
    ## if the no. of columns are 3, expect "sample", "barcode", "cluster", if 4, expect sample:bc, umap1, umap2, cluster
    df = pd.read_csv(args.groupInfo, sep="\t", index_col=None, comment="#")
    if len(df.columns) == 3:
        df.columns = ["sample", "barcode", "cluster"]
    elif len(df.columns) == 4:
        df.columns = ["Cell_ID", "umap1", "umap2", "cluster"]
        try:
            df[["sample", "barcode"]] = df.Cell_ID.str.split("::", expand=True)
        except:
            sys.exit(df_err)
    else:
        sys.exit(df_err)

    df.index = df[["sample", "barcode"]].apply(lambda x: "::".join(x), axis=1)
    # barcodes = groupInfo.barcode.unique().tolist()

    ## match the sample labels with groupInfo labels before proceeding
    ## create new DF with user-provided bam+labels (this would match the counts)

    # check that sample names in df and --labels match
    labels_not_in_df = set(args.labels).difference(set(df["sample"]))
    df_not_in_labels = set(df["sample"]).difference(set(args.labels))
    if bool(df_not_in_labels):
        sys.stderr.write(
            "Some (or all) of the samples indicated in the groupInfo file "
            "are absent from the bam file labels! \n"
            "Mismatched samples are: {} \n".format(df_not_in_labels)
        )
        df = df.loc[df["sample"].isin(set(args.labels)),]
    elif bool(labels_not_in_df):
        sys.stderr.write(
            "Some (or all) of the samples indicated in --labels "
            "are absent from the in the groupInfo file! \n"
            "Mismatched samples are: {} \n".format(labels_not_in_df)
        )

    # make another df, containing union of all barcodes, and in the same order per sample
    barcodes = df["barcode"].unique().tolist()
    sm = list(itertools.chain.from_iterable(itertools.repeat(lab, len(barcodes)) for lab in args.labels))
    bc = barcodes * len(args.labels)
    groupInfo = pd.DataFrame({"sample": sm, "barcode": bc})
    groupInfo.index = groupInfo[["sample", "barcode"]].apply(lambda x: "::".join(x), axis=1)
    groupInfo = pd.merge(
        groupInfo,
        df["cluster"],
        how="left",
        left_index=True,
        right_index=True,
        sort=False,
    )
    groupInfo = groupInfo.reset_index()[["sample", "barcode", "cluster"]]
    groupInfo["cluster"] = groupInfo["cluster"].astype("category")
    # re-construct new labels (sample+bc)
    newlabels = ["::".join([x, y]) for x, y in zip(groupInfo["sample"], groupInfo["barcode"])]
    # Normalization options
    scale_factor = args.scaleFactor
    func_args = {"scaleFactor": args.scaleFactor}

    coverageAsFrequency = False
    if not args.ignoreForNormalization:
        args.ignoreForNormalization = []
    if args.normalizeUsing == "None":
        args.normalizeUsing = None  # For the sake of sanity
    elif args.normalizeUsing == "Frequency":
        coverageAsFrequency = True
        args.normalizeUsing = "Mean"  # mean of binarized counts shall give us frequency

    # This fixes issue #520, where --extendReads wasn't honored if --filterRNAstrand was used
    if args.filterRNAstrand and not args.Offset:
        args.Offset = [1, -1]

    if args.MNase:
        # check that library is paired end
        # using getFragmentAndReadSize
        from deeptools.getFragmentAndReadSize import get_read_and_fragment_length

        fraglengths = [
            get_read_and_fragment_length(
                file,
                return_lengths=False,
                blackListFileName=args.blackListFileName,
                numberOfProcessors=args.numberOfProcessors,
                verbose=args.verbose,
            )
            for file in args.bamfiles
        ]
        if any([x[0] is None for x in fraglengths]):
            sys.exit("*Error*: For the --MNAse function a paired end library is required. ")

        # Set some default fragment length bounds
        if args.minFragmentLength == 0:
            args.minFragmentLength = 130
        if args.maxFragmentLength == 0:
            args.maxFragmentLength = 200

        wr = CenterFragment(
            args.bamfiles,
            binLength=args.binSize,
            stepSize=args.binSize,
            barcodes=barcodes,
            clusterInfo=groupInfo,
            cellTag=args.cellTag,
            groupTag=args.groupTag,
            groupLabels=newlabels,
            motifFilter=args.motifFilter,
            genome2bit=args.genome2bit,
            GCcontentFilter=args.GCcontentFilter,
            region=args.region,
            blackListFileName=args.blackListFileName,
            numberOfProcessors=args.numberOfProcessors,
            extendReads=args.extendReads,
            minMappingQuality=args.minMappingQuality,
            duplicateFilter=args.duplicateFilter,
            center_read=args.centerReads,
            zerosToNans=False,
            samFlag_include=args.samFlagInclude,
            samFlag_exclude=args.samFlagExclude,
            minFragmentLength=args.minFragmentLength,
            maxFragmentLength=args.maxFragmentLength,
            chrsToSkip=args.ignoreForNormalization,
            binarizeCoverage=coverageAsFrequency,
            verbose=args.verbose,
        )

    elif args.Offset:
        if len(args.Offset) > 1:
            if args.Offset[0] == 0:
                sys.exit(
                    "*Error*: An offset of 0 isn't allowed, since offsets are 1-based positions inside each alignment."
                )
            if args.Offset[1] > 0 and args.Offset[1] < args.Offset[0]:
                sys.exit("'Error*: The right side bound is less than the left-side bound. This is inappropriate.")
        else:
            if args.Offset[0] == 0:
                sys.exit(
                    "*Error*: An offset of 0 isn't allowed, since offsets are 1-based positions inside each alignment."
                )
        wr = OffsetFragment(
            args.bamfiles,
            binLength=args.binSize,
            stepSize=args.binSize,
            barcodes=barcodes,
            clusterInfo=groupInfo,
            cellTag=args.cellTag,
            groupTag=args.groupTag,
            groupLabels=newlabels,
            motifFilter=args.motifFilter,
            genome2bit=args.genome2bit,
            GCcontentFilter=args.GCcontentFilter,
            region=args.region,
            numberOfProcessors=args.numberOfProcessors,
            extendReads=args.extendReads,
            minMappingQuality=args.minMappingQuality,
            duplicateFilter=args.duplicateFilter,
            center_read=args.centerReads,
            zerosToNans=False,
            samFlag_include=args.samFlagInclude,
            samFlag_exclude=args.samFlagExclude,
            minFragmentLength=args.minFragmentLength,
            maxFragmentLength=args.maxFragmentLength,
            chrsToSkip=args.ignoreForNormalization,
            binarizeCoverage=coverageAsFrequency,
            verbose=args.verbose,
        )
        wr.filter_strand = args.filterRNAstrand
        wr.Offset = args.Offset
    else:
        wr = WriteBedGraph.WriteBedGraph(
            args.bamfiles,
            binLength=args.binSize,
            stepSize=args.binSize,
            barcodes=barcodes,
            clusterInfo=groupInfo,
            cellTag=args.cellTag,
            groupTag=args.groupTag,
            groupLabels=newlabels,
            motifFilter=args.motifFilter,
            genome2bit=args.genome2bit,
            GCcontentFilter=args.GCcontentFilter,
            region=args.region,
            blackListFileName=args.blackListFileName,
            numberOfProcessors=args.numberOfProcessors,
            extendReads=args.extendReads,
            minMappingQuality=args.minMappingQuality,
            duplicateFilter=args.duplicateFilter,
            center_read=args.centerReads,
            zerosToNans=False,
            samFlag_include=args.samFlagInclude,
            samFlag_exclude=args.samFlagExclude,
            minFragmentLength=args.minFragmentLength,
            maxFragmentLength=args.maxFragmentLength,
            chrsToSkip=args.ignoreForNormalization,
            binarizeCoverage=coverageAsFrequency,
            verbose=args.verbose,
        )

    wr.run(
        WriteBedGraph.scaleCoverage,
        func_args,
        args.outFilePrefix,
        blackListFileName=args.blackListFileName,
        normUsing=args.normalizeUsing,
        format=args.outFileFormat,
        smoothLength=None,
    )


class OffsetFragment(WriteBedGraph.WriteBedGraph):
    """
    Class to redefine the get_fragment_from_read for the --Offset case
    """

    def filterStrand(self, read, rv):
        """
        A generic read filtering function that gets used by everything in this class.

        rv is returned if the strand is correct, otherwise [(None, None)]
        """
        # Filter by RNA strand, if desired
        if read.is_paired:
            if self.filter_strand == "forward":
                if read.flag & 144 == 128 or read.flag & 96 == 64:
                    return rv
            elif self.filter_strand == "reverse":
                if read.flag & 144 == 144 or read.flag & 96 == 96:
                    return rv
            else:
                return rv
        else:
            if self.filter_strand == "forward":
                if read.flag & 16 == 16:
                    return rv
            elif self.filter_strand == "reverse":
                if read.flag & 16 == 0:
                    return rv
            else:
                return rv

        return [(None, None)]

    def get_fragment_from_read_list(self, read, offset):
        """
        Return the range of exons from the 0th through 1st bases, inclusive. Positions are 1-based
        """
        rv = [(None, None)]
        blocks = read.get_blocks()
        blockLen = sum([x[1] - x[0] for x in blocks])

        if self.defaultFragmentLength != "read length":
            if self.is_proper_pair(read, self.maxPairedFragmentLength):
                if read.is_reverse:
                    foo = (read.next_reference_start, read.reference_start)
                    if foo[0] < foo[1]:
                        blocks.insert(0, foo)
                else:
                    foo = (
                        read.reference_end,
                        read.reference_end + abs(read.template_length) - read.infer_query_length(),
                    )
                    if foo[0] < foo[1]:
                        blocks.append(foo)

            # Extend using the default fragment length
            else:
                if read.is_reverse:
                    foo = (
                        read.reference_start - self.defaultFragmentLength + read.infer_query_length(),
                        read.reference_start,
                    )
                    if foo[0] < 0:
                        foo = (0, foo[1])
                    if foo[0] < foo[1]:
                        blocks.insert(0, foo)
                else:
                    foo = (
                        read.reference_end,
                        read.reference_end + self.defaultFragmentLength - read.infer_query_length(),
                    )
                    if foo[0] < foo[1]:
                        blocks.append(foo)

        stretch = []
        # For the sake of simplicity, convert [(10, 20), (30, 40)] to [10, 11, 12, 13, ..., 40]
        # Then subset accordingly
        for block in blocks:
            stretch.extend(range(block[0], block[1]))
        if read.is_reverse:
            stretch = stretch[::-1]

        # Handle --centerReads
        if self.center_read:
            _ = (len(stretch) - blockLen) // 2
            stretch = stretch[_ : _ + blockLen]

        # Subset by --Offset
        try:
            foo = stretch[offset[0] : offset[1]]
        except:
            return rv

        if len(foo) == 0:
            return rv
        if read.is_reverse:
            foo = foo[::-1]

        # Convert the stretch back to a list of tuples
        foo = np.array(foo)
        d = foo[1:] - foo[:-1]
        idx = np.argwhere(d > 1).flatten().tolist()  # This now holds the interval bounds as a list
        idx.append(-1)
        last = 0
        rv = []
        for i in idx:
            rv.append((foo[last].astype("int"), foo[i].astype("int") + 1))
            last = i + 1

        # Handle strand filtering, if needed
        return self.filterStrand(read, rv)

    def get_fragment_from_read(self, read):
        """
        This is mostly a wrapper for self.get_fragment_from_read_list(),
        which needs a list and for the offsets to be tweaked by 1.
        """
        offset = [x for x in self.Offset]
        if len(offset) > 1:
            if offset[0] > 0:
                offset[0] -= 1
            if offset[1] < 0:
                offset[1] += 1
        else:
            if offset[0] > 0:
                offset[0] -= 1
                offset = [offset[0], offset[0] + 1]
            else:
                if offset[0] < -1:
                    offset = [offset[0], offset[0] + 1]
                else:
                    offset = [offset[0], None]
        if offset[1] == 0:
            # -1 gets switched to 0, which screws things up
            offset = (offset[0], None)
        return self.get_fragment_from_read_list(read, offset)


class CenterFragment(WriteBedGraph.WriteBedGraph):
    """
    Class to redefine the get_fragment_from_read for the --MNase case

    The coverage of the fragment is defined as the 2 or 3 basepairs at the
    center of the fragment length.
    """

    def get_fragment_from_read(self, read):
        """
        Takes a proper pair fragment of high quality and limited
        to a certain length and outputs the center
        """
        fragment_start = fragment_end = None

        # only paired forward reads are considered
        # Fragments have already been filtered according to length
        if read.is_proper_pair and not read.is_reverse and 1 < abs(read.tlen):
            # distance between pairs is even return two bases at the center
            if read.tlen % 2 == 0:
                fragment_start = read.pos + read.tlen / 2 - 1
                fragment_end = fragment_start + 2

            # distance is odd return three bases at the center
            else:
                fragment_start = read.pos + read.tlen / 2 - 1
                fragment_end = fragment_start + 3

        return [(fragment_start, fragment_end)]
