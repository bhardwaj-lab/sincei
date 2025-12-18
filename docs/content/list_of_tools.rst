Command Line tools
==================

.. toctree::
    :maxdepth: 1
    :hidden:

    tools/scFilterBarcodes
    tools/scFilterStats
    tools/scJSD
    tools/scCountReads
    tools/scCountQC
    tools/scCombineCounts
    tools/scClusterCells
    tools/scBulkCoverage
    tools/scFindVCRs

The following tools use **BAM files** as input. These BAM files could can from any single-cell genomics protocol, as long as they have a **tag** that specifies the cell barcodes.

* :doc:`tools/scFilterBarcodes`
* :doc:`tools/scFilterStats`
* :doc:`tools/scCountReads`
* :doc:`tools/scBulkCoverage`
* :doc:`tools/scJSD`

The following tools use the `AnnData <https://anndata.readthedocs.io/>`_ output produced within the sincei analysis workflow, with file extension **.h5ad**.

* :doc:`tools/scCountQC`
* :doc:`tools/scCombineCounts`
* :doc:`tools/scClusterCells`
* :doc:`tools/scFindVCRs`

The table below summarizes input and output of each tool:

.. contents::
    :local:

+--------------------------------------+------------------+-----------------------+---------------------------------------------+-------------------------------------------------------------------------------------------------------------+
| tool                                 | type             | input files           | main output file(s)                         | application                                                                                                 |
+======================================+==================+=======================+=============================================+=============================================================================================================+
|:ref:`scFilterBarcodes`               | preprocessing    | BAM/SAM files         | text file with filtered cell barcodes       | Identify and filter cell barcodes (for droplet-based single-cell seq)                                       |
+--------------------------------------+------------------+-----------------------+---------------------------------------------+-------------------------------------------------------------------------------------------------------------+
|:ref:`scFilterStats`                  | QC               | BAM/SAM files         | text file with QC per cell                  | Produce per-cell statistics after filtering reads by user-defined criteria.                                 |
+--------------------------------------+------------------+-----------------------+---------------------------------------------+-------------------------------------------------------------------------------------------------------------+
|:ref:`scCountReads`                   | preprocessing    | BAM/SAM files         | h5ad object with cellxregion counts         | Counts reads for each barcode on genomic bins or user-defined features.                                     |
+--------------------------------------+------------------+-----------------------+---------------------------------------------+-------------------------------------------------------------------------------------------------------------+
|:ref:`scCountQC`                      | QC               | h5ad object           | QC metrics / filtered h5ad file             | Perform quality control and filter the output of scCountReads.                                              |
+--------------------------------------+------------------+-----------------------+---------------------------------------------+-------------------------------------------------------------------------------------------------------------+
|:ref:`scCombineCounts`                | preprocessing    | h5ad objects          | merged h5ad object                          | Concatenate/merge the counts from different samples/batches or modalities                                   |
+--------------------------------------+------------------+-----------------------+---------------------------------------------+-------------------------------------------------------------------------------------------------------------+
|:ref:`scClusterCells`                 | analysis         | h5ad object           | .tsv file with clusters, png with UMAP      | Perform dimensionality reduction and clustering on the output of scCountReads.                              |
+--------------------------------------+------------------+-----------------------+---------------------------------------------+-------------------------------------------------------------------------------------------------------------+
|:ref:`scBulkCoverage`                 | analysis         | tsv file + BAM file   | bigwig files                                | Get pseudo-bulk coverage per group using a user-supplied cell->group mapping (output of scClusterCells).    |
+--------------------------------------+------------------+-----------------------+---------------------------------------------+-------------------------------------------------------------------------------------------------------------+
|:ref:`scFindVCRs`                     | analysis         | h5ad object           | bed files                                   | Call variable chromatin regions (VCRs) from binned chromatin data.                                          |
+--------------------------------------+------------------+-----------------------+---------------------------------------------+-------------------------------------------------------------------------------------------------------------+