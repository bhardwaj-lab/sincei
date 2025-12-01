.. image:: ./content/images/sincei-logo-transparent.png
    :align: left
    :width: 350

.. image:: https://readthedocs.org/projects/sincei/badge/?version=latest
    :target: https://sincei.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://img.shields.io/pypi/v/sincei.svg?style=plastic
    :target: https://pypi.org/project/sincei/
    :alt: PyPI Version

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
    :target: http://bioconda.github.io/recipes/sincei/README.html
    :alt: Install with bioconda


=======================
Single Cell Informatics
=======================

A user-friendly toolkit for QC, counting, clustering and plotting of single-cell (epi)genomics data.

sincei is described in our preprint: `Bhardwaj V. , Mourragui, S. (2024) User-friendly exploration of epigenomic data in single cells using sincei. <https://www.biorxiv.org/content/10.1101/2024.07.27.605424v1>`__


============
Installation
============

sincei is a command line toolkit based on python3. The stable version of sincei can be installed using `conda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ , while the development versions can be installed from github via pip.

Installation via bioconda
-------------------------

Create a new conda environment and install sincei using:

.. code-block:: bash

    conda create -n sincei -c bioconda -c conda-forge sincei

Mac users with Arm architecture (M1-3 macbooks and beyond) should explicitly specify the osx-64 version to allow dependencies to install properly.

.. code-block:: bash

    conda create -n sincei --subdir 'osx-64' -c bioconda -c conda-forge sincei

*Note:* The dependency `mctorch-lib` required for `scClusterCells` is currently unavailable via conda, therefore, to use `scClusterCells`, we recommend installing it separately via pip.

.. code-block:: bash

    # install mctorch-lib
    (sincei): pip install mctorch-lib
    (sincei): scClusterCells --help


Installation via github
-----------------------

Create a new conda environment and install sincei in it using pip from GitHub:

.. code-block:: bash

    conda create -n sincei python=3.10
    conda activate sincei
    (sincei): pip install git+https://github.com/bhardwaj-lab/sincei.git@master#egg=sincei


Getting Help
------------

* For questions related to usage, or suggesting changes/enhancements please use our `GitHub discussion board <https://github.com/bhardwaj-lab/sincei/discussions>`__ . To report bugs, please create an issue on `our GitHub repository <https://github.com/bhardwaj-lab/sincei>`_

**Please Note that sincei is under active development.** We expect significant changes/updates as we move towards our first major release (1.0).


Command line tools available in sincei
--------------------------------------

Tools for a typical single-cell analysis workflow (WIP: work in progress/not available yet)

========================== ============================================================================================================
tool                                 description
========================== ============================================================================================================
:ref:`scFilterBarcodes`        Identify and filter cell barcodes from BAM file (for droplet-based single-cell seq).
:ref:`scFilterStats`           Produce per-cell statistics after filtering reads by user-defined criteria.
:ref:`scJSD`                   Compare the cumulative read coverages of cells on regions using the Jensen-Shannon distance.
:ref:`scCountReads`            Counts reads for each barcode on genomic bins or user-defined features.
:ref:`scCountQC`               Perform quality control and filter the output of scCountReads.
:ref:`scCombineCounts`         Concatenate/merge the counts from different samples/batches or different modalities.
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
