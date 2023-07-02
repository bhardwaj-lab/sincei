from itertools import compress
from deeptools.utilities import getTLen
import numpy as np
import sys


def checkBAMtag(bam, name, tag):
    bctag = [read.has_tag(tag) for read in bam.head(1000)]
    if not any(bctag):
        sys.stderr.write("WARNING: Input file {} seems to lack the tag {}. Output might be empty. \n".format(name, tag))
    return None


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
    read_motif = read.get_forward_sequence()[0 : len(readMotif)]
    ref_motifLen = len(refMotif) - 1

    if read.is_reverse:
        # for reverse reads ref motif begins at read-end and ends downstream
        endpos = read.reference_end + ref_motifLen
        if endpos > genome.chroms()[chrom]:
            endpos = read.reference_end  # fail without error
        ref_motif = genome.sequence(chrom, read.reference_end - 1, endpos)
    else:
        # for forward reads ref motif begins upstream and ends at read-start
        startpos = read.reference_start - ref_motifLen
        if startpos < 0:
            startpos = read.reference_start  # fail without error
        ref_motif = genome.sequence(chrom, startpos, read.reference_start + 1)

    if read_motif == readMotif and ref_motif == refMotif:
        return True
    else:
        return False


def checkGCcontent(read, lowFilter, highFilter, returnGC=False):
    r"""Checks if the GC content of the read is within the given range

    Parameters
    ----------
    read : pysam.AlignedSegment
        A pysam AlignedSegment object

    lowFilter : float
        Minimum GC content

    highFilter : float
        Maximum GC content

    returnGC : bool
        If true, return the GC content of the read

    Returns
    -------
    bool
        True if the GC content of the read is within the given range

    Examples
    --------

    >>> test = Tester()
    >>> read = test.bamFile1.fetch().next()
    >>> checkGCcontent(read, 0.3, 0.7)
    True
    """
    seq = read.get_forward_sequence()
    total_bases = len(seq)
    gc_bases = len([x for x in seq if x == "C" or x == "G"])
    gc_frac = float(gc_bases) / total_bases
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
    matchSum = sum([i[1] for i in list(compress(cig, matchPos))])
    if matchSum / tot >= lowFilter:
        return True
    else:
        return False


def colorPicker(name):
    r"""
    This function returns a list of colors for plotting.

    Parameters
    ----------
    name : str
        The name of the color palette to use.

    Returns
    -------
    list
        A list of colors.

    Examples
    --------

    >>> colorPicker('twentyfive')
    ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999']

    >>> colorPicker('colorblind')
    ['#0072B2', '#009E73', '#D55E00', '#CC79A7', '#F0E442', '#56B4E9', '#E69F00', '#000000']
    """

    colors = {
        "twentyfive": [
            "#dd4e34",
            "#78d545",
            "#b047dc",
            "#d6d94d",
            "#5d48c6",
            "#74d88b",
            "#c944af",
            "#8a9739",
            "#542c75",
            "#d3953c",
            "#607dc7",
            "#487f46",
            "#d04774",
            "#7bd2c8",
            "#6b2737",
            "#cfcb9a",
            "#332a42",
            "#d7928a",
            "#343d25",
            "#cc8ad3",
            "#7b6a43",
            "#b5bad9",
            "#99472b",
            "#4e8290",
            "#936987",
        ],
        "colorblind": [
            "#ffaec0",
            "#cd7600",
            "#893558",
            "#195f37",
            "#da71f9",
            "#a2d391",
            "#881e9b",
            "#b9d05f",
            "#524d7e",
            "#f2bf4b",
            "#01d6c2",
            "#f54040",
            "#0097fb",
            "#756400",
            "#4b44aa",
            "#f0bd79",
            "#a1008b",
            "#6d4d02",
            "#ff9afc",
            "#01df8f",
            "#e2b8ed",
            "#6e9d00",
            "#f4177b",
            "#01b65f",
            "#9b2532",
        ],
    }

    return colors[name]


def getDupFilterTuple(read, bc, filterArg):
    r"""
    Returns a tuple with the information needed to filter duplicates, based on read and filter type.
    The tuple is composed of the barcode, the umi, the start and end positions
    and the chromosome name.

    Parameters
    ----------
    read : pysam.AlignedSegment
        A pysam.AlignedSegment object

    bc : str
        The barcode

    filter : str
        A string with the type of filter to use.

    Returns
    -------
    tuple

        A tuple with the information needed to filter duplicates.
        The tuple is composed of the barcode, the umi, the start and end positions
        and the chromosome name.

    Examples
    --------

    >>> test = Tester()
    >>> read = test.bamFile1.fetch().next()
    >>> getDupFilterTuple(read, 'ATCG', 'end_umi')
    ('ATCG', 'ATCG', None, None, 0, False)
    """

    tLenDup = getTLen(read, notAbs=True)
    filt = filterArg.split("_")

    # get fragment start and end for that read
    if tLenDup >= 0:
        s = read.pos
        e = s + tLenDup
    else:
        s = read.pnext
        e = s - tLenDup

    if "end" not in filt:
        # ignore read/fragment end and mate information
        mate_refid = read.reference_id
        if read.is_reverse:
            s = None
        else:
            e = None
    else:
        # use mate info, reset fragment end to mate pos if read is chimeric
        if read.reference_id != read.next_reference_id:
            e = read.pnext
        mate_refid = read.next_reference_id

    # get UMI if asked
    if "umi" in filt:
        umi = read.get_tag("RX")
    else:
        umi = None

    tup = (bc, umi, s, e, mate_refid, read.is_reverse)
    return tup


def gini(i, X):
    r"""Computes the Gini coefficient for each row of a sparse matrix (Obs*Var).

    Parameters
    ----------
    i : int
        row index
    X : numpy array
        matrix

    Returns
    -------
    float
        Gini coefficient for the given row

    Examples
    --------

    >>> X = np.matrix([[1,2,3,4],[5,6,7,8],[9,10,11,12]])
    >>> gini(0, X)
    0.0
    >>> gini(1, X)
    0.0
    >>> gini(2, X)
    0.0
    """

    # based on bottom eq:
    # http://www.statsdirect.com/help/generatedimages/equations/equation154.svg
    # from:
    # http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    # All values are treated equally, arrays must be 1d:
    array = X[i, :].A.flatten()  # get all bins from i'th cell
    array = array[array.nonzero()]

    if array.shape[0] <= 1:
        return np.nan
    else:
        array = np.sort(array)
        # Index per array element:
        index = np.arange(1, array.shape[0] + 1)
        # Number of array elements:
        n = array.shape[0]
        # Gini coefficient:
        return (np.sum((2 * index - n - 1) * array)) / (n * np.sum(array))
