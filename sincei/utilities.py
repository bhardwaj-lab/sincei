from itertools import compress
from deeptools.utilities import getTLen

def checkMotifs(read, chrom, genome, readMotif, refMotif):
        """
        Check whether a given motif is present in the read and the corresponding reference genome.
        For example, in MNAse (scChIC-seq) data, we expect the reads to have an 'A' at the 5'-end,
        while the genome has a 'TA' over hang (where the 'A' aligns with 'A' in the forward read),
        like this below.

        Forwards aligned read: read has 'A', upstream has T
        R1 ........A------->
        ----------TA------------\ Ref (+)

        Rev aligned read: read has 'T', downstream has A

        <-------T....... R1
        --------TA------------\ Ref (+)

        This function can look for any arbitrary motif in read and corresponding genome, but in the
        same orientation as described above.

        :return: bool

        >>> import pysam
        >>> import os
        >>> from scDeepTools.scReadCounter import CountReadsPerBin as cr
        >>> root = os.path.dirname(os.path.abspath(__file__)) + "/test/test_data/"
        >>> bam = pysam.AlignmentFile("{}/test_TA_filtering.bam".format(root))
        >>> iter = bam.fetch()
        >>> read = next(iter)
        >>> cr.checkMotifs(read, 'A', 'TA') # for valid scChIC read
        True
        >>> read = next(iter)
        >>> cr.checkMotifs(read, 'A', 'TA') # for invalid scChIC read
        False
        """
        # get read and ref motif pos
        read_motif = read.get_forward_sequence()[0:len(readMotif)]
        ref_motifLen = len(refMotif) - 1

        if(read.is_reverse):
            # for reverse reads ref motif begins at read-end and ends downstream
            endpos = read.reference_end + ref_motifLen
            if endpos > genome.chroms()[chrom]:
                endpos = read.reference_end #fail without error
            ref_motif = genome.sequence(chrom, read.reference_end - 1, endpos)
        else:
            # for forward reads ref motif begins upstream and ends at read-start
            startpos = read.reference_start - ref_motifLen
            if startpos < 0:
                startpos = read.reference_start #fail without error
            ref_motif = genome.sequence(chrom, startpos, read.reference_start + 1)

        if read_motif == readMotif and ref_motif == refMotif:
            return True
        else:
            return False


def checkGCcontent(read, lowFilter, highFilter, returnGC=False):
    """
    Check whether the read falls into the range of min and max GC content provided.
    Return: Bool

    """
    seq = read.get_forward_sequence()
    total_bases = len(seq)
    gc_bases = len([x for x in seq if x == 'C' or x == 'G'])
    gc_frac = float(gc_bases)/total_bases
    if returnGC:
        return gc_frac
    else:
        if gc_frac >= lowFilter and gc_frac <= highFilter:
            return True
        else:
            return False


def checkAlignedFraction(read, lowFilter):
    """
    Check whether the fraction of read length that aligns to the reference is higher than
    the given threshold. Aligned fraction includes the max allowed mismatches tolerated by
    the aligner, and excludes InDels and Clippings.

    Return: Bool
    """
    cig = read.cigartuples
    tot = read.infer_read_length()
    matchPos = [i[0] == 0 for i in cig]
    matchSum = sum([i[1] for i in list(compress(cig, matchPos)) ])
    if matchSum/tot >= lowFilter:
        return True
    else:
        return False

def colorPicker(name):
    """
    select colors from my list of custom made color color_palettes
    """
    colors = {'twentyfive': ["#dd4e34",  "#78d545",  "#b047dc",  "#d6d94d",  "#5d48c6",  "#74d88b",  "#c944af",
                 "#8a9739",  "#542c75",  "#d3953c",  "#607dc7",  "#487f46",  "#d04774",  "#7bd2c8",
                 "#6b2737",  "#cfcb9a",  "#332a42",  "#d7928a",  "#343d25",  "#cc8ad3",  "#7b6a43",
                 "#b5bad9",  "#99472b",  "#4e8290",  "#936987"],
              'colorblind': ["#ffaec0", "#cd7600", "#893558", "#195f37", "#da71f9", "#a2d391", "#881e9b",
                 "#b9d05f", "#524d7e", "#f2bf4b", "#01d6c2", "#f54040", "#0097fb", "#756400",
                 "#4b44aa", "#f0bd79", "#a1008b", "#6d4d02", "#ff9afc", "#01df8f", "#e2b8ed",
                 "#6e9d00", "#f4177b", "#01b65f", "#9b2532"]
            }

    return colors[name]


def getDupFilterTuple(read, bc, filter):
    """
    based on read and filter type, return a tuple to match with previous read, such that duplicates can be detected.
    """
    tLenDup = getTLen(read, notAbs=True)
    filt = filter.split('_')
    ## get read (or fragment) start/end
    # get fragment start and end for that read
    if tLenDup >= 0:
        s = read.pos
        e = s + tLenDup
    else:
        s = read.pnext
        e = s - tLenDup
    if read.reference_id != read.next_reference_id:
        e = read.pnext
    if 'end' not in filt:
        # use only read (or fragment) start
        if read.is_reverse:
            s = None
        else:
            e = None
    ## get UMI if asked
    if 'umi' in filt:
        umi = read.get_tag('RX')
    else:
        umi = None
    tup = (bc, umi, s, e, read.next_reference_id, read.is_reverse)
    return tup
