The CLI tools
==================

.. toctree::
    :maxdepth: 1
    :hidden:

    tools/scFilterBarcodes
    tools/scFilterStats
    tools/scCountReads
    tools/scCountQC
    tools/scCombineCounts
    tools/scClusterCells
    tools/scBulkCoverage
    tools/scFindVCRs
    tools/scScoreFeatures
    tools/scPlotRegion
    tools/scExportSignal

The following tools use **BAM files** as input. These BAM files could can from any
single-cell genomics protocol, as long as they have a **tag** that specifies the cell
barcodes:

* :doc:`tools/scFilterBarcodes`
* :doc:`tools/scFilterStats`
* :doc:`tools/scCountReads`
* :doc:`tools/scBulkCoverage`

The following tools use the `AnnData <https://anndata.readthedocs.io/>`_ output produced
by other sincei commands, with file extension **.h5ad**:

* :doc:`tools/scCountQC`
* :doc:`tools/scCombineCounts`
* :doc:`tools/scClusterCells`
* :doc:`tools/scFindVCRs`
* :doc:`tools/scScoreFeatures`
* :doc:`tools/scPlotRegion`
* :doc:`tools/scExportSignal`

The table below summarizes input and output of each tool:

.. contents::
    :local:

.. list-table::
   :header-rows: 1
   :widths: 18 12 15 22 33

   * - tool
     - type
     - input files
     - main output file(s)
     - application
   * - :ref:`scFilterBarcodes`
     - preprocessing
     - :ref:`BAM files <fmt-bam>`
     - :ref:`text file <fmt-tsv>` with filtered cell barcodes
     - Identify and filter cell barcodes (for droplet-based single-cell seq).
   * - :ref:`scFilterStats`
     - QC
     - :ref:`BAM files <fmt-bam>`
     - :ref:`text file <fmt-tsv>` with QC per cell
     - Produce per-cell statistics after filtering reads by user-defined criteria.
   * - :ref:`scCountReads`
     - preprocessing
     - :ref:`BAM files <fmt-bam>`
     - :ref:`h5ad object <fmt-anndata>` with cellxregion counts
     - Counts reads for each barcode on genomic bins or user-defined features.
   * - :ref:`scCountQC`
     - QC
     - :ref:`h5ad object <fmt-anndata>`
     - QC metrics (:ref:`tsv <fmt-tsv>`) / filtered :ref:`h5ad file <fmt-anndata>`
     - Perform quality control and filter the output of scCountReads.
   * - :ref:`scFindVCRs`
     - analysis
     - :ref:`h5ad object <fmt-anndata>`
     - :ref:`BED file <fmt-bed>`
     - Call variable chromatin regions (VCRs) from binned chromatin data.
   * - :ref:`scScoreFeatures`
     - analysis
     - :ref:`h5ad object <fmt-anndata>` + :ref:`BED <fmt-bed>`/:ref:`GTF <fmt-gtf>`
     - :ref:`h5ad object <fmt-anndata>` with gene activity scores
     - Calculate gene activity scores from chromatin features/bins or aggregate bins
       into VCRs.
   * - :ref:`scCombineCounts`
     - preprocessing
     - :ref:`h5ad objects <fmt-anndata>`
     - merged :ref:`h5ad <fmt-anndata>` or :ref:`h5mu <fmt-mudata>` object
     - Concatenate/merge the counts from different samples/batches or modalities.
   * - :ref:`scClusterCells`
     - analysis
     - :ref:`h5ad object <fmt-anndata>`
     - :ref:`.tsv file <fmt-tsv>` with clusters, :ref:`png <fmt-plots>` with UMAP
     - Perform dimensionality reduction and clustering on the output of scCountReads.
   * - :ref:`scBulkCoverage`
     - analysis
     - :ref:`tsv file <fmt-tsv>` + :ref:`BAM file <fmt-bam>`
     - :ref:`bigwig files <fmt-bigwig>`
     - Get pseudo-bulk coverage per group using a user-supplied cell->group mapping
       (output of scClusterCells).
   * - :ref:`scPlotRegion`
     - analysis
     - :ref:`h5ad object <fmt-anndata>`
     - :ref:`png <fmt-plots>`
     - Plot genomic regions of interest from sincei output.
   * - :ref:`scExportSignal`
     - analysis
     - :ref:`h5ad object <fmt-anndata>`
     - :ref:`bedgraph matrix <fmt-bm>`, :ref:`mtx <fmt-mtx>`, loom
     - Export sincei-supportred .h5ad (anndata) object to other formats.
