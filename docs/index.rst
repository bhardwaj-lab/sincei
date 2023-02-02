
sincei
=======

.. image:: ./content/images/sincei.png

Bhardwaj V. (2022) sincei: A user-friendly toolkit for QC, counting, clustering and plotting of single-cell (epi)genomics data.


Installation
-------------------
sincei is a command line toolkit based on python3, and can be installed using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Create a new conda environment and install sincei using:

.. code-block:: bash

    cd <programs_folder>
    conda create -n sincei -c bioconda -c conda-forge scanpy deeptools
    conda activate sincei
    (sincei): pip install --editable=git+https://github.com/vivekbhr/sincei.git@master#egg=sincei



The list of tools available in sincei
---------------------------------------

Tools for a typical single-cell analysis workflow (WIP: work in progress/not available yet)

========================== ============================================================================================================
tool                                 description
========================== ============================================================================================================
:ref:`scFilterBarcodes`        Identify and filter cell barcodes from BAM file (for droplet-based single-cell seq)
:ref:`scFilterStats`           Produce per-cell statistics after filtering reads by user-defined criteria.
:ref:`scCountReads`            Counts reads for each barcode on genomic bins or user-defined features.
:ref:`scCountQC`               Perform quality control and filter the output of scCountReads.
:ref:`scCombineCounts`         Concatenate/merge the counts from different samples/batches or modalities
:ref:`scClusterCells`          Perform dimensionality reduction and clustering on the output of scCountReads.
:ref:`scBulkCoverage`          Get pseudo-bulk coverage per group using a user-supplied cell->group mapping (output of scClusterCells).
      scFindMarkers            [WIP] Find marker genes per group, given the output of scCountReads and a user-defined group.
      scFeaturePlot            [WIP] Plot the counts for a given feature on a UMAP or on a (IGV-style) genomic-track.
========================== ============================================================================================================


Getting Help
------------

* For all kind of questions, suggesting changes/enhancements and to report bugs, please create an issue on `our GitHub repository <https://github.com/deeptools/HiCExplorer>`_


Contents:
---------

.. toctree::
   :maxdepth: 2

   content/list_of_tools
   content/tutorials
   content/modules/modules
   content/news


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
