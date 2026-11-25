|pypi| |bioconda| |test| |docs| |codecov| |license| |black|

.. |pypi| image:: https://img.shields.io/pypi/v/sincei.svg?style=plastic
   :alt: PyPI Version
   :target: https://pypi.org/project/sincei/
.. |bioconda| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
   :alt: install with bioconda
   :target: http://bioconda.github.io/recipes/sincei/README.html
.. |test| image:: https://github.com/vivekbhr/sincei/actions/workflows/test.yml/badge.svg
   :alt: test
.. |docs| image:: https://readthedocs.org/projects/sincei/badge/?version=latest
   :alt: Documentation Status
   :target: https://sincei.readthedocs.io/en/latest/?badge=latest
.. |codecov| image:: https://codecov.io/gh/vivekbhr/sincei/graph/badge.svg?token=VRTMITHHBI
   :alt: codecov
   :target: https://codecov.io/gh/vivekbhr/sincei
.. |license| image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :alt: License: MIT
   :target: https://opensource.org/licenses/MIT
.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :alt: Code style: black
   :target: https://github.com/psf/black

.. image:: ./content/images/sincei-logo-light.png
    :align: left
    :class: only-light main-logo

.. image:: ./content/images/sincei-logo-dark.png
    :align: left
    :class: only-dark main-logo

====================================================================
A command-line toolkit for exploring single-cell (epi)genomics data.
====================================================================

Sincei is a flexible command-line toolkit designed to work with single-cell data directly
from BAM files:

- Perform read-level and count-level quality control, and produce
  `MultiQC <https://github.com/multiqc/multiqc>`_ compatible reports.
- Make count and coverage matrices from single-cell BAM files, for any single-cell genomics
  protocol.
- Aggregate signal in bins, genes, peaks, or any feature of interest from single-cells.
- Perform dimensionality reduction and clustering of all kinds of single-cell data (open
  chromatin, histone marks, methylation, gene expression, etc.).
- Create coverage files (bigwigs) per group of cells for visualization.
- Along with additional tools to visualize, interpret, and annotate cells and genomic regions.

sincei is also part of the `scverse <https://scverse.org/>`_ ecosystem, and it's command-line tools can easily be combined with other Python or R packages for further analysis.

For details, please `read our preprint <https://www.biorxiv.org/content/10.1101/2024.07.27.605424v1>`_ describing sincei.


Getting Help
------------

For questions related to usage, or suggesting changes and new features please use our
`GitHub discussion board <https://github.com/bhardwaj-lab/sincei/discussions>`__ .

To report bugs, please create an issue on
`our GitHub repository <https://github.com/bhardwaj-lab/sincei>`_

**Please Note that sincei is under active development.** We expect significant changes/updates
as we move towards our first major release (1.0).


Command line tools available in sincei
--------------------------------------

Tools for a typical single-cell analysis workflow

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - tool
     - description
   * - :ref:`scFilterBarcodes`
     - Identify and filter cell barcodes from BAM file (for droplet-based single-cell seq).
   * - :ref:`scFilterStats`
     - Produce per-cell statistics after filtering reads by user-defined criteria.
   * - :ref:`scCountReads`
     - Counts reads for each barcode on genomic bins or user-defined features.
   * - :ref:`scCountQC`
     - Perform quality control and filter the output of scCountReads.
   * - :ref:`scFindVCRs`
     - Call variable chromatin regions (VCRs) from binned chromatin data.
   * - :ref:`scScoreFeatures`
     - Calculate gene activity scores from chromatin features/bins.
   * - :ref:`scCombineCounts`
     - Concatenate/merge the counts from different samples/batches or different modalities.
   * - :ref:`scClusterCells`
     - Perform dimensionality reduction and clustering on the output of scCountReads.
   * - :ref:`scBulkCoverage`
     - Get pseudo-bulk coverage per group using a user-supplied cell->group mapping (output of scClusterCells).
   * - :ref:`scPlotRegion`
     - Plot genomic regions of interest from sincei output.
   * - :ref:`scExportSignal`
     - Export sincei-supportred .h5ad (anndata) object to other formats.


Contents
---------

.. toctree::
    :maxdepth: 2

    content/installation
    content/list_of_tools
    content/tutorials
    content/modules
    content/methods
    content/file_formats
    content/references
    content/news
