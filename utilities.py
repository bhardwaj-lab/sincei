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
        >>> cr.is_proper_pair(read, 'A', 'TA') # for invalid scChIC read
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


def checkGCcontent(read, lowFilter, highFilter):
    """
    Check whether the read falls into the range of min and max GC content provided.
    Return: Bool
    
    """
    seq = read.get_forward_sequence()
    total_bases = len(seq)
    gc_bases = len([x for x in seq if x == 'C' or x == 'G'])
    gc_frac = float(gc_bases)/total_bases
    if gc_frac >= lowFilter and gc_frac <= highFilter:
        return True
    else:
        return False
