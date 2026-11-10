File formats
============

This page describes the file formats that sincei reads and writes.

Input files
-----------

.. _fmt-bam:

BAM (.bam, and index .bai)
~~~~~~~~~~~~~~~~~~~~~~~~~~

BAM is the compressed, binary form of SAM (Sequence Alignment Map), the standard format
to store aligned sequencing reads. It is the input of the sincei tools that process
sequencing reads. A BAM file has a header, followed by one record per alignment. View it
as text with ``samtools view``.

The **header** (lines that start with ``@``) lists the chromosomes and their lengths
(``@SQ``) and the read groups (``@RG``, optional). The :ref:`scCountReads` option
``--groupTag`` requires ``@RG`` lines.

Each **alignment record** has 11 required fields and optional tags, separated by
tabs. sincei uses the following fields:

.. list-table::
   :header-rows: 1
   :widths: 10 25 65

   * - Field
     - Name
     - Use in sincei
   * - 2
     - FLAG
     - Bits that describe properties of the alignment, for example 16 (reverse
       strand), 64 (first mate), 256 (secondary) or 2048 (supplementary). Here is a
       `useful tool <https://broadinstitute.github.io/picard/explain-flags.html>`_
       to check the meaning of SAM flags combinations.
   * - 3
     - RNAME
     - The chromosome name of the alignment.
   * - 4
     - POS
     - The 1-based start position of the alignment
   * - 5
     - MAPQ
     - The mapping quality, the confidence of the aligner in the position of the read.
   * - 6
     - CIGAR
     - The matches, insertions, deletions and clipped bases of the alignment.
   * - 9
     - TLEN
     - The template length, the length of a read pair's fragment of origin.
   * - 12+
     - Tags
     - Optional fields such as with a two-letter name, a type and a
       value. sincei requires input BAMs to have a cell barcode tag and, optionally,
       a UMI tag, a sample tag (for merged BAMs with many samples) or a tag with a
       value to count per read.

Input BAM files must be **sorted by position** and **indexed**. The index is a ``.bai``
file with the same name and in the same directory as the BAM file (``sample.bam.bai``),
which lets sincei read regions in the BAM without needing to read the whole file. Make
it with ``samtools index sample.bam``. sincei only reads BAI indexes, not CSI indexes.

.. note::

   The exact definition of the mapping quality, the CIGAR string and the optional tags
   depend on the aligner used to produce the BAM. Check how your aligner defines them
   before you use them to filter reads with sincei.

.. _fmt-bed:

BED (.bed)
~~~~~~~~~~

BED (Browser Extensible Data) is a text format for genomic regions, with one region
per line and tab-separated columns:

.. code-block:: text

   chr1    3000000    3010000    region_1
   chr1    3500000    3502000    region_2

- Columns 1 to 3 (required): the chromosome, the start and the end. The start is
  0-based and the end is not included (half-open interval notation), so the first
  region above covers 10,000 bases.
- Column 4 (optional): the name of the region. sincei uses it in the region IDs.
  If there is no name, sincei uses ``chrom:start-end``.
- More columns are allowed, consult
  `this <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_ for more details.

BED lines that start with ``#`` are skipped. The file can be compressed with gzip or
bgzip. sincei reads BED files of features to count (:ref:`scCountReads`,
:ref:`scScoreFeatures`) and blacklists (``--blacklist``).

:ref:`scFindVCRs` writes a 6-column BED file (``chrom``, ``start``, ``end``,
``name``, ``score``, ``strand``). The score column is the value of the penalty
parameter used for the segmentation.

.. _fmt-gtf:

GTF and GFF3 (.gtf, .gff3)
~~~~~~~~~~~~~~~~~~~~~~~~~~

GTF (Gene Transfer Format) and GFF3 (Generic Feature Format version 3) are text formats
for gene and genome annotation, with nine tab-separated columns. Unlike BED, the
coordinates are 1-based and the end is included.

.. code-block:: text

   chr1  HAVANA  gene  3073253  3074322  .  +  .  gene_id "ENSMUSG00000102693"; gene_name "4933401J01Rik";

sincei uses these columns:

- Column 1: the chromosome. Columns 4 and 5: the start and the end. Column 7: the
  strand.
- Column 3: the feature type. ``--featureID`` (default ``gene``) selects the features
  to count, and ``--exonID`` (default ``exon``) their exons.
- Column 9: the attributes. ``--featureIDtag`` selects the attribute that names
  each feature: by default ``gene_id`` for GTF and ``ID`` for GFF3.

The file can be compressed with gzip or bgzip. Download annotations from, for
example, GENCODE or Ensembl, for the same genome assembly as the BAM files.

Consult the
`GTF <https://biocorecrg.github.io/PhD_course_genomics_format_2021/gtf_format.html>`_
and
`GFF3 <https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md>`_
specifications for more details.

.. _fmt-2bit:

2bit (.2bit)
~~~~~~~~~~~~

2bit is a compressed binary format for a genome sequence. ``--motifFilter`` needs it
(``--genome2bit``) to read the genome overhang of reads. Download it from UCSC, or
make it from a FASTA file with ``faToTwoBit genome.fa genome.2bit``.

.. _fmt-barcodes:

Cell barcode lists (.txt)
~~~~~~~~~~~~~~~~~~~~~~~~~

A cell barcode list is a plain text file with one barcode per line:

.. code-block:: text

   AAACAGCCAAGGAATC
   AAACAGCCAATCCCTT

``--barcodes`` gives the barcodes to keep (a whitelist). :ref:`scFilterStats` needs
one. :ref:`scFilterBarcodes` reads one with ``--whitelist``. The cell blacklist of
:ref:`scCountQC` (``--cellBlacklist``) is a similar list, with one cell ID
(``sample::barcode``) per line.

Output files
------------

.. _fmt-anndata:

AnnData (.h5ad)
~~~~~~~~~~~~~~~

`AnnData <https://anndata.scverse.org/en/stable/>`_ is the data structure that sincei
uses to store count matrices, saved as ``.h5ad`` files (based on HDF5). Read it in
Python with ``anndata.read_h5ad("counts.h5ad")`` or with
`scanpy <https://scanpy.readthedocs.io/en/stable/index.html>`_.
It is compatible with other `scverse <http://scverse.org/>`_ tools and many
`bioconductor <https://bioconductor.org/>`_ tools via
`anndataR <https://anndatar.scverse.org/index.html>`_.

An AnnData object has the following parts:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Part
     - Content in sincei
   * - ``X``
     - The count matrix: one row per cell and one column per feature. It may be a dense
       or a sparse matrix. Sincei always stores data in sparse format.
   * - ``obs``, ``obs_names``
     - A table with one row per cell. The cell IDs are ``<sample>::<barcode>``.
       sincei needs the columns ``sample`` and ``barcode``. :ref:`scCountQC` adds
       per-cell metrics, and :ref:`scClusterCells` adds clusters (``leiden``).
   * - ``var``, ``var_names``
     - A table with one row per region. Sincei names features as
       ``<chrom>_<start>_<end>::<name>`` (``::None`` for bins). sincei needs the
       columns ``chrom``, ``start``, ``end``, and ``name``. :ref:`scCountQC` adds
       per-region metrics.
   * - ``obsm``
     - Matrices with one row per cell, for example the dimensionality reductions of
       :ref:`scClusterCells` (``X_glmpca``, ``X_pca``, ``X_lsa``, ``X_lda``) and
       UMAP coordinates (``X_umap``).
   * - ``varm``
     - Matrices with one row per region, for example the loadings of
       :ref:`scClusterCells` (``glmpca``, ``pca``, ``lsa``, ``lda``).
   * - ``uns``
     - Other information, for example the parameters of a method.

``--compression`` (default ``none``, or ``gzip``) of :ref:`scCountReads` can compress
the file with gzip.

.. _fmt-mudata:

MuData (.h5mu)
~~~~~~~~~~~~~~

MuData is the Python data structure for multi-modal data, saved as ``.h5mu`` files.
It holds one AnnData object for each modality (for example chromatin and RNA) of the
same cells. :ref:`scCombineCounts` writes it with ``--method multi-modal``. Read it
with ``mudata.read_h5mu`` or the `muon <https://muon.scverse.org/stable/#>`_ package.
See the `MuData documentation <https://mudata.readthedocs.io/>`_.

.. _fmt-tsv:

Tables (.tsv)
~~~~~~~~~~~~~

Several tools write tab-separated tables with a header line:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Tool
     - Columns
   * - :ref:`scFilterBarcodes`
     - ``barcode``, ``count`` (the number of non-zero bins), ``selected``,
       ``count_log10``, ``count_rank``.
   * - :ref:`scFilterStats`
     - ``Cell_ID``, ``Total_sampled`` (a count), and the percentage of the sampled
       reads that each filter removes (for example ``Low_MAPQ`` or
       ``Internal_Duplicates``). MultiQC can read it.
   * - :ref:`scCountQC`
     - ``<prefix>.cells.tsv``: ``Cell_ID`` and the cell metrics.
       ``<prefix>.regions.tsv``: ``Feature_ID`` and the region metrics.
   * - :ref:`scClusterCells`
     - ``Cell_ID``, ``UMAP1``, ``UMAP2``, ``cluster``. Written next to the UMAP plot.

:ref:`scBulkCoverage` reads a table of cell groups (``--groupInfo``): either the
4-column output of :ref:`scClusterCells`, or a 3-column table with ``sample``,
``barcode`` and ``group``.

.. _fmt-bigwig:

bigWig and bedGraph (.bw, .bedgraph)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:ref:`scBulkCoverage` writes the pseudo-bulk coverage of each cell group as
``<prefix>_<group>.bw`` (default) or ``<prefix>_<group>.bedgraph``
(``--outFileFormat``).

- **bedGraph** is a text format with four columns: chromosome, start, end (as in
  BED) and a value, for example the coverage.
- **bigWig** is the binary, indexed form of bedGraph. It is small and fast to read,
  so use it in genome browsers such as IGV or the UCSC Genome Browser. Consult its
  `specification <https://genome.ucsc.edu/goldenpath/help/bigWig.html>`_ for more
  details.

.. _fmt-bm:

BedGraphMatrix (.bm)
~~~~~~~~~~~~~~~~~~~~

BedGraphMatrix (``.bm``) is a plain text format with one row per region and one column
per cell, without a header. The first three columns are the chromosome, the start and
the end, and the other columns are the counts of each cell. The rows are sorted by
position. :ref:`scExportSignal` writes it (``--outFileFormat bm``) for
pyGenomeTracks, which plots it as a heatmap of cells.

.. _fmt-mtx:

MatrixMarket (.mtx)
~~~~~~~~~~~~~~~~~~~

MatrixMarket (``.mtx``) is a palin text format to store sparse matrices: it lists only
the non-zero entries, with their row and column. :ref:`scExportSignal` writes three
files:

- ``<prefix>.counts.mtx``: the counts, with one row per cell and one column per
  region.
- ``<prefix>.rownames.txt``: the cell IDs, one per line.
- ``<prefix>.colnames.txt``: the region IDs, one per line.

R (``Matrix::readMM``) and many single-cell tools can read it.

.. _fmt-plots:

Plots (.png, .jpg, .svg, .pdf)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The tools that make plots write them as ``png`` (default), ``jpg``, ``svg`` or
``pdf`` (``--plotFileFormat``), with ``--plotWidth`` and ``--plotHeight`` in cm and
``--dpi`` (default 300).
