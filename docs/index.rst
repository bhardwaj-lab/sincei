

.. image:: ./content/images/sincei-logo.png
      :align: left
      :width: 350

.. image:: https://zenodo.org/badge/271841139.svg
    :target: https://zenodo.org/badge/latestdoi/271841139
    :alt: DOI

.. image:: https://readthedocs.org/projects/sincei/badge/?version=latest)]
    :target: https://sincei.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://img.shields.io/pypi/v/sincei.svg?style=plastic
    :target: https://pypi.org/project/sincei/
    :alt: PyPI Version

Bhardwaj V. and Mourragui S. (2023) sincei: A user-friendly toolkit for QC, counting, clustering and plotting of single-cell (epi)genomics data.




Installation
-------------------
sincei is a command line toolkit based on python3, and can be installed using `conda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_

Create a new conda environment and install sincei using:

.. code-block:: bash

    conda create -n sincei -c anaconda python=3.8
    conda activate sincei
    (sincei): pip install git+https://github.com/vivekbhr/sincei.git@master#egg=sincei


Getting Help
------------

* For all kind of questions, suggesting changes/enhancements or to report bugs, please create an issue on `our GitHub repository <https://github.com/vivekbhr/sincei>`_

**Please Note that sincei is under active development.** Some features might be incomplete, untested or might be removed as we move towards a stable version.



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


Contents:
---------

.. toctree::
   :maxdepth: 2

   content/list_of_tools
   content/tutorials
   content/modules/modules
   content/news
