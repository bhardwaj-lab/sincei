Creating a count matrix
=======================

:ref:`scCountReads` makes a count matrix :math:`X \in \mathbb{R}^{N \times P}` from one
or more BAM files, with one row per cell and one column per feature. The features are
either bins from a uniform tiling of the genome, or regions from an annotation file.
The entry :math:`x_{ci}` is the number of reads of cell :math:`c` that overlap feature
:math:`i`.

Every read from a BAM file goes through the same steps in both modes:

1. The filters decide if a read is counted or filtered out.
2. The read is turned into a genomic interval (multiple intervals if spliced).
3. The interval is adjusted (e.g., read extension, centering, Tn5 correction).
4. The interval adds to the features it overlaps.

We must note an of step 4, it can produce either:

- A *coverage matrix*, where reads are counted in every feature they overlap, and thus
  a read may contribute more than one value.
- A *count matrix*, where reads are counted only once, in the feature with the largest
  overlap.

By default, sincei produces coverage matrices. Section `Coverage and count matrices`_
explains how to produce each of them.

Filters
-------

sincei applies the filters to each read in this order. A read that fails one filter
is not counted.

1. **Flags and mapping quality.** Unmapped reads are removed. Optional filters remove
   reads with a low mapping quality, with required or disallowed SAM flags (multiple
   SAM flags can be specified at the same time), or select RNA-seq reads by the strand
   of their gene.
2. **Region.** When a region is specified, only reads in this region are counted.
3. **Blacklist.** A read is removed when blacklisted regions cover at least half of
   its alignment span. The filter uses the alignment span of the whole read, not the
   interval of step 2.
4. **Barcode.** A read without the barcode tag is removed. When a barcode whitelist is
   given, reads with barcodes not in the list are filtered out.
5. **Read properties.** Optional limits on fragment length, GC content and aligned
   fraction of each read. The fragment length is the template length
   :math:`|\mathrm{TLEN}|` of the BAM record for paired reads, and the reference bases
   that the read covers for single-end reads.
6. **Duplicates.** Optional, see below.
7. **Motif.** Optional. Only reads with a motif at its 5' end and reference genome
   overhang (for example the cut site of a restriction enzyme) are kept. Multiple 
   motifs can be specified.

Duplicates
~~~~~~~~~~

The duplicate filter keeps the first read of each group of duplicates. Two reads are
considered duplicates based on one of the following criteria:

- ``start_bc``: reads with the same start and cell barcode are considered duplicates.
- ``start_bc_umi``: reads with the same start, cell barcode, and UMI are considered
  duplicates.
- ``start_end_bc``: reads with the same start, end, and cell barcode are considered
  duplicates.
- ``start_end_bc_umi``: reads with the same start, end, cell barcode, and UMI are
  considered duplicates.

``start_bc_umi`` and ``start_end_bc_umi`` require unique molecular identifiers (UMIs)
to be present in the input BAM files in a tag (``--umiTag``).

The filter always compares reads that start at the same base. Since input BAM files are
sorted and processed in chunks, we keep a small set of duplicate keys in memory at a
time.

For paired-end data, the two mates of a fragment are different reads are both kept.

From a read to an interval
--------------------------

By default, a read adds its alignment span :math:`[a, e)` (0-based, end not included).
A spliced read adds its aligned blocks, so introns get no counts. The following options
modify the interval:

- ``--extendReads``: reads are extended to the length of their fragment. For a properly
  paired read that faces its mate, on the same chromosome, and with a fragment length
  not more than four times the extension length, the fragment goes from the beginning of
  a read to its the end of its mate (TLEN). Other reads are extended in their direction
  to a fixed length. Without a specified extension length, sincei uses the mean of the
  median fragment lengths of the BAM files, from their first 100,000 reads.
- ``--centerReads``: the interval becomes a window as wide as the read is at the center
  of the fragment.
- ``--offset``: only some aligned bases of the read count, for example the first base
  (1) or the last base (-1). Two values can be specified to take a part of the read
  (e.g, ``--offset 1 5`` takes the 5 bases a read).
- ``--tn5``: perform Tn5 correction on the data, the fragment becomes
  :math:`[a + 4, e - 5)`, to remove the 9-base duplication that Tn5 makes. Only the
  forward mate of a proper pair is counted, so each fragment counts once.
- ``--mnase``: to determine nucleosome positions from MNase data, only the two center
  bases of each fragment count. By default, the minimum and maximum fragment lengths
  are set to 130 to 200 bases, to remove dinucleosomes. Only the forward mate of a
  proper pair is counted.

``--tn5`` and ``--mnase`` centering require paired-end data.

Combining ``--centerReads`` with ``--offset`` gives a window of the size specified by
``--offset`` at the center of the fragment.

Coverage and count matrices
---------------------------

Sincei produces coverage matrices by default, where a read is counted in every feature
it overlaps. When counting on a non-overlapping genome tiling or feature set, we can
use ``--offset`` flag to count each read only once. For example:

- ``--offset 1`` counts reads where their first base is aligned.
- ``--offset -1`` counts reads where their last base is aligned.
- ``--offset 1 --centerReads`` counts reads in the center of the read.
- ``--offset 1 --extendReads --centerReads --samFlagInclude 64`` counts fragments
  where their center is. We only keep the first read of a fragment
  (``--samFlagInclude 64``), extend it to its fragment size (``--extendReads``), and
  count it where its middle base falls (``--offset 1 --centerReads``).

Counting into bins
------------------

The bin counting mode tiles each chromosome of the first input BAM file. With bin size
:math:`b` and the step :math:`s` between the starts of two bins, bin :math:`j` of a
chromosome of length :math:`L` is:

.. math::

   [\,s j,\; \min(s j + b,\, L)\,), \qquad j = 0, \dots, \lceil L / s \rceil - 1

By default :math:`s = b`, so the bins touch but do not overlap, and the last bin of a
chromosome may be shorter. With :math:`s > b` there are gaps between the bins. When a
region is given, the bins tile only that region only.

A read counts once in every bin that its interval overlaps by default. A read that
crosses the border between two bins adds one to each of them. Thus, the sum of a row
can be larger than the number of reads of the cell. The shorter the bins, the more
often this will happen.

Unlike tools like `deeptools bamCoverage
<https://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html>`_,
sincei does not exclude the regions in a blacklist (``--blacklist``) from the tiling.
Deeptools bins skip over the blacklisted regions, meaning not all the starts of a
region will be multiples of the specified bin size. Sincei guarantees that the start of
every bin remains a multiple of the bin size.

Counting into features
----------------------

The features counting mode count a feature set from a BED, GTF or GFF3 file (see
:doc:`../file_formats`). A read counts once in every feature that its interval overlaps,
also when the features overlap each other.

For scRNA-seq data, we recommend using the ``--metagene`` flag. It counts only in the
exons of a gene or transcript. A read then counts only once: the gene/transcript with
the most aligned bases in common with the read gets the count. Thus, a read that
overlaps two exons of one gene counts once, and a read in the intron of a gene does not
count. We recommend using a GTF or GFF3 with the gene annotations of a reference genome
when using ``--metagene``.

When a region is given, the features that start in the region are kept, also when
they continue after its end.

Values and cells
----------------

By default, each read adds 1 to the features it is counted in. Alternatively, a read
can add the value of a numeric tag in its BAM (``--valueTag``). Reads without the tag
are not counted.

Cell IDs
--------

The rows in a count matrix are cells, with their ID set to ``sample::barcode``. The
sample is the name of the BAM file. For merged BAMs with a group tag (``--groupTag``)
one BAM file can hold several samples. The header of the merged BAM file must declare
the list of sample names.

When a barcode whitelist is given, the matrix will hold one row for each sample and
barcode in the list, regardless of whether the "cell" has no counts. Without a list,
the matrix has one row for each barcode with at least one count in a BAM file.

Bear in mind a row may not contain the counts of a real cell, simply those of a barcode
in the dataset, specially for droplet data. Pre-filtering barcodes from a whitelist
using :ref:`scFilterBarcodes` and filtering the matrix produced by :ref:`scCountReads`
with :ref:`scCountQC` allows us to identify both real cells and high quality cells for
analysis.

Parallel counting
-----------------

sincei divides each chromosome into chunks (1 Mb by default), produces a list of chunks
per bam, and counts each chunk as a job. A BAM is only read if one of its chunks is
currently being processed, capping the memory used at a time.

Each chunk reads its part of the BAM file through the BAI index. A read belongs to the
chunk that contains its alignment start, so no read is counted twice. The interval of a
read can continue into the next chunk, so the counts of two chunks for the same cell
and feature are added to an accumulator.

To make efficient use of memory, the results are added to the output AnnData file when
all the chunks of a BAM are finished.

Sparsity and storage
--------------------

Single-cell data is inherently sparse, often with >90% of zero entries in a matrix.
Sincei stores the results of :ref:`scCountReads` in
`compressed sparse row <http://en.wikipedia.org/wiki/Sparse_matrix>`_ (CSR) format
to greatly reduce the disk space needed to store the data matrix.

The results of :ref:`scCountReads` can be stored in compressed format using
``--compression gzip``. This also helps in saving disk space at the cost of having to
decompress data before use.

Space and time complexity
-------------------------

With :math:`R` total reads in the BAM files:

- **Time**: :math:`O(R)` for bin counting, because the bins of a read come from its
  position. For feature counting, :math:`O(R \log P)`, because each read searches an
  interval index of the features on its chromosome.
- **Memory**: :math:`O(\mathrm{nnz})` for the data matrix, plus the counts of the
  chunks in progress. The output is written one BAM file at a time.
