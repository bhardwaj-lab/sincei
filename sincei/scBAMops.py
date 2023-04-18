#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pysam
import os
import sys
import pandas as pd
import py2bit
from deeptools import parserCommon
from deeptools.bamHandler import openBam
from deeptools.mapReduce import mapReduce
from deeptools.utilities import getTLen, smartLabels, getTempFileName
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

# own functions
# scriptdir=os.path.abspath(os.path.join(__file__, "../../sincei"))
# sys.path.append(scriptdir)

from sincei import ParserCommon
from sincei.Utilities import checkMotifs, checkGCcontent, getDupFilterTuple
from sincei._version import __version__

## UPDATE: add group tag to BAM file based on a 2-columns mapping file (barcode -> group)


def parseArguments():
    internalParser = get_args()
    ioParser = ParserCommon.inputOutputOptions(opts=["bamfile", "groupInfo", "outFile"])
    bamParser = ParserCommon.bamOptions(suppress_args=["binSize", "distanceBetweenBins"])
    filterParser = ParserCommon.filterOptions()
    readParser = ParserCommon.readOptions(suppress_args=["extendReads", "centerReads"])
    otherParser = ParserCommon.otherOptions()
    parser = argparse.ArgumentParser(
        parents=[
            ioParser,
            internalParser,
            bamParser,
            filterParser,
            readParser,
            otherParser,
        ],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="This tool filters alignments in a BAM/CRAM file according the the specified parameters. It can optionally output to BEDPE format.",
        usage="Example usage: scBAMops -b sample1.bam -o sample1.filtered.bam --minMappingQuality 10 --filterMetrics log.txt",
        add_help=False,
    )
    return parser


def get_args():
    parser = argparse.ArgumentParser(add_help=False)

    mod = parser.add_argument_group("Read modifiction arguments")
    mod.add_argument(
        "--shift",
        nargs="+",
        type=int,
        help='Shift the left and right end of a read (for BAM files) or a fragment (for BED files). A positive value shift an end to the right (on the + strand) and a negative value shifts a fragment to the left. Either 2 or 4 integers can be provided. For example, "2 -3" will shift the left-most fragment end two bases to the right and the right-most end 3 bases to the left. If 4 integers are provided, then the first and last two refer to fragments whose read 1 is on the left or right, respectively. Consequently, it is possible to take strand into consideration for strand-specific protocols. A fragment whose length falls below 1 due to shifting will not be written to the output. See the online documentation for graphical examples. Note that non-properly-paired reads will be filtered.',
    )

    mod.add_argument(
        "--ATACshift",
        action="store_true",
        help="Shift the produced BAM file or BEDPE regions as commonly done for ATAC-seq. This is equivalent to --shift 4 -5 5 -4.",
    )

    output = parser.add_argument_group("Optional output arguments")
    output.add_argument(
        "--BED",
        action="store_true",
        help="Instead of producing BAM files, write output in BEDPE format (as defined by MACS2). Note that only reads/fragments passing filtering criterion are written in BEDPE format.",
    )

    output.add_argument(
        "--filterMetrics",
        metavar="FILE.log",
        help="The number of entries in total and filtered are saved to this file",
    )

    output.add_argument(
        "--filteredOutReads",
        metavar="filtered.bam",
        help="If desired, all reads NOT passing the filtering criteria can be written to this file.",
    )

    output.add_argument(
        "--outTagName",
        metavar="STR",
        type=str,
        help="In case where you want to group the BAM files based on user-defined --groupInfo, specify the output BAM tag which would contain the group name.",
        default=None,
    )

    return parser


def shiftRead(b, chromDict, args):
    if not b.is_proper_pair:
        return None
    tLen = getTLen(b, notAbs=True)
    start = b.pos
    end = start + b.query_alignment_end
    if b.is_reverse and not b.is_read2:
        end -= args.shift[2]
        deltaTLen = args.shift[3] - args.shift[2]
    elif b.is_reverse and b.is_read2:
        end += args.shift[1]
        deltaTLen = args.shift[1] - args.shift[0]
    elif not b.is_reverse and not b.is_read2:
        start += args.shift[0]
        deltaTLen = args.shift[1] - args.shift[0]
    else:
        start -= args.shift[3]
        deltaTLen = args.shift[3] - args.shift[2]

    # Sanity check
    if end - start < 1:
        if b.is_reverse:
            start = end - 1
        else:
            end = start + 1
    if start < 0:
        start = 0
    if end > chromDict[b.reference_name]:
        end = chromDict[b.reference_name]
    if end - start < 1:
        return None

    # create a new read
    b2 = pysam.AlignedSegment()
    b2.query_name = b.query_name
    b2.flag = b.flag
    b2.reference_id = b.reference_id
    b2.reference_start = start
    b2.mapping_quality = b.mapping_quality
    b2.cigar = ((0, end - start),)  # Returned cigar is only matches
    if tLen < 0:
        b2.template_length = tLen - deltaTLen
    else:
        b2.template_length = tLen + deltaTLen
    b2.next_reference_id = b.next_reference_id
    b2.next_reference_start = b.next_reference_start
    if b.is_proper_pair:
        if b2.is_read2 and b2.is_reverse:
            b2.next_reference_start += args.shift[0]
        elif not b2.is_read2 and b2.is_reverse:
            b2.next_reference_start -= args.shift[3]

    return b2


def filterWorker(arglist):
    chrom, start, end, args, chromDict = arglist
    fh = openBam(args.bamfile)
    # open 2 bit if needed
    if args.genome2bit:
        twoBitGenome = py2bit.open(args.genome2bit, True)
        if chrom not in twoBitGenome.chroms().keys():
            raise NameError("chromosome {} not found in 2bit file".format(chrom))

    mode = "wbu"
    oname = getTempFileName(suffix=".bam")
    ofh = pysam.AlignmentFile(oname, mode=mode, template=fh)

    if args.filteredOutReads:
        onameFiltered = getTempFileName(suffix=".bam")
        ofiltered = pysam.AlignmentFile(onameFiltered, mode=mode, template=fh)
    else:
        onameFiltered = None
        ofiltered = None

    prev_pos = set()
    lpos = None

    nFiltered = 0
    total = 0
    for read in fh.fetch(chrom, start, end):
        if read.pos < start:
            # ensure that we never double count (in case distanceBetweenBins == 0)
            continue
        total += 1

        bc = read.get_tag(args.cellTag)
        if isinstance(args.groupInfo, pd.DataFrame):
            smpl = read.get_tag(args.groupTag)
            try:
                tag_to_add = args.groupInfo.loc[
                    (args.groupInfo["sample"] == smpl) & (args.groupInfo["barcode"] == bc),
                    "cluster",
                ].values.item()
                read.set_tag(args.outTagName, tag_to_add, value_type="Z", replace=True)
            except:
                if args.verbose:
                    sys.stderr.write("Encountered read tags not in groupInfo file, skipped..")
                nFiltered += 1
                if ofiltered:
                    ofiltered.write(read)
                continue

        ## -- begin filtering -- ##
        if read.flag & 4:
            # Ignore unmapped reads, they were counted already
            nFiltered += 1
            if ofiltered:
                ofiltered.write(read)
            continue

        if args.minMappingQuality and read.mapq < args.minMappingQuality:
            nFiltered += 1
            if ofiltered:
                ofiltered.write(read)
            continue

        if args.samFlagInclude and read.flag & args.samFlagInclude != args.samFlagInclude:
            nFiltered += 1
            if ofiltered:
                ofiltered.write(read)
            continue
        if args.samFlagExclude and read.flag & args.samFlagExclude != 0:
            nFiltered += 1
            if ofiltered:
                ofiltered.write(read)
            continue

        tLen = getTLen(read)
        if args.minFragmentLength > 0 and tLen < args.minFragmentLength:
            nFiltered += 1
            if ofiltered:
                ofiltered.write(read)
            continue
        if args.maxFragmentLength > 0 and tLen > args.maxFragmentLength:
            nFiltered += 1
            if ofiltered:
                ofiltered.write(read)
            continue

        if args.duplicateFilter:
            tup = getDupFilterTuple(read, bc, args.duplicateFilter)
            if lpos is not None and lpos == read.reference_start and tup in prev_pos:
                nFiltered += 1
                if ofiltered:
                    ofiltered.write(read)
                continue
            if lpos != read.reference_start:
                prev_pos.clear()
            lpos = read.reference_start
            prev_pos.add(tup)

        # remove reads with low/high GC content
        if args.GCcontentFilter:
            if not checkGCcontent(read, args.GCcontentFilter[0], args.GCcontentFilter[1]):
                nFiltered += 1
                if ofiltered:
                    ofiltered.write(read)
                continue

        # remove reads that don't pass the motif filter
        if args.motifFilter:
            if not checkMotifs(read, chrom, twoBitGenome, args.motifFilter[0], args.motifFilter[1]):
                nFiltered += 1
                if ofiltered:
                    ofiltered.write(read)
                continue

        # filterRNAstrand
        if args.filterRNAstrand:
            if read.is_paired:
                if args.filterRNAstrand == "forward":
                    if read.flag & 144 == 128 or read.flag & 96 == 64:
                        pass
                    else:
                        nFiltered += 1
                        if ofiltered:
                            ofiltered.write(read)
                        continue
                elif args.filterRNAstrand == "reverse":
                    if read.flag & 144 == 144 or read.flag & 96 == 96:
                        pass
                    else:
                        nFiltered += 1
                        if ofiltered:
                            ofiltered.write(read)
                        continue
            else:
                if args.filterRNAstrand == "forward":
                    if read.flag & 16 == 16:
                        pass
                    else:
                        nFiltered += 1
                        if ofiltered:
                            ofiltered.write(read)
                        continue
                elif args.filterRNAstrand == "reverse":
                    if read.flag & 16 == 0:
                        pass
                    else:
                        nFiltered += 1
                        if ofiltered:
                            ofiltered.write(read)
                        continue

        if args.shift:
            read = shiftRead(read, chromDict, args)
            if not read:
                continue

        # Read survived filtering
        ofh.write(read)

    # The results from the workers will get sorted, so get the TID
    tid = fh.get_tid(chrom)

    ofh.close()
    if ofiltered:
        ofiltered.close()
    fh.close()
    if args.genome2bit:
        twoBitGenome.close()

    return tid, start, total, nFiltered, oname, onameFiltered


def convertBED(oname, tmpFiles, chromDict):
    """
    Stores results in BEDPE format, which is:
    chromosome	frag_leftend	frag_rightend

    The fragment ends can be shifted
    """
    ofile = open(oname, "w")
    for tmpFile in tmpFiles:
        fh = pysam.AlignmentFile(tmpFile)

        for b in fh.fetch(until_eof=True):
            tLen = getTLen(b, notAbs=True)
            if tLen > 0:
                start = b.pos
                end = start + tLen
                if end > chromDict[b.reference_name]:
                    end = chromDict[b.reference_name]
                if end - start < 1:
                    continue
                ofile.write("{}\t{}\t{}\n".format(b.reference_name, start, end))
        fh.close()
        os.unlink(tmpFile)
    ofile.close()


def main(args=None):
    args = parseArguments().parse_args(args)

    # grouping and tagging the output alignments
    if args.groupInfo:
        if not args.outTagName:
            sys.stderr.write(
                "--groupInfo provided without --outTagName. By default, groups will be added to the 'RG' tag in the output BAM file."
            )
            args.outTagName = "RG"

        df = pd.read_csv(args.groupInfo, sep="\t", index_col=None, comment="#")
        if len(df.columns) == 3:
            df.columns = ["sample", "barcode", "cluster"]
            df.index = df[["sample", "barcode"]].apply(lambda x: "::".join(x), axis=1)
            args.groupInfo = df
        else:
            sys.stderr.write(
                "Error: The input to --groupInfo must be a 3-column tsv file with <sample> <barcode> and <group>"
            )

    # shifting the output alignments
    if args.shift:
        if len(args.shift) not in [2, 4]:
            sys.exit("The --shift option can accept either 2 or 4 values only.")
        if len(args.shift) == 2:
            args.shift.extend([-args.shift[1], -args.shift[0]])
    elif args.ATACshift:
        args.shift = [4, -5, 5, -4]

    bam, mapped, unmapped, stats = openBam(args.bamfile, returnStats=True, nThreads=args.numberOfProcessors)
    total = mapped + unmapped
    chrom_sizes = [(x, y) for x, y in zip(bam.references, bam.lengths)]
    chromDict = {x: y for x, y in zip(bam.references, bam.lengths)}

    # if motifFilter is given, convert the input to list
    if args.motifFilter:
        if not args.genome2bit:
            print("MotifFilter asked but genome (2bit) file not provided.")
            sys.exit(1)
        else:
            args.motifFilter = args.motifFilter.strip(" ").split(",")

    if args.GCcontentFilter:
        gc = args.GCcontentFilter.strip(" ").split(",")
        args.GCcontentFilter = [float(x) for x in gc]

    # Filter, writing the results to a bunch of temporary files
    res = mapReduce(
        [args, chromDict],
        filterWorker,
        chrom_sizes,
        region=args.region,
        blackListFileName=args.blackListFileName,
        numberOfProcessors=args.numberOfProcessors,
        verbose=args.verbose,
    )

    res = sorted(res)  # The temp files are now in order for concatenation
    nFiltered = sum([x[3] for x in res])
    totalSeen = sum([x[2] for x in res])  # The * contig isn't queried

    tmpFiles = [x[4] for x in res]
    if not args.BED:
        arguments = ["-o", args.outFile]
        arguments.extend(tmpFiles)  # [..., *someList] isn't available in python 2.7
        pysam.samtools.cat(*arguments)
        for tmpFile in tmpFiles:
            os.unlink(tmpFile)
    else:
        convertBED(args.outFile, tmpFiles, chromDict)

    if args.filteredOutReads:
        tmpFiles = [x[5] for x in res]
        if not args.BED:
            arguments = ["-o", args.filteredOutReads]
            arguments.extend(tmpFiles)  # [..., *someList] isn't available in python 2.7
            pysam.samtools.cat(*arguments)
            for tmpFile in tmpFiles:
                os.unlink(tmpFile)
        else:
            convertBED(args.outFile, tmpFiles, chromDict, args)

    if args.filterMetrics:
        sampleName = args.bamfile
        if args.label:
            sampleName = args.label
        if args.smartLabels:
            sampleName = smartLabels([args.bamfile])[0]

        of = open(args.filterMetrics, "w")
        of.write("#bamFilterReads --filterMetrics\n")
        of.write("#File\tReads Remaining\tTotal Initial Reads\n")
        of.write("{}\t{}\t{}\n".format(sampleName, totalSeen - nFiltered, total))
        of.close()

    return 0
