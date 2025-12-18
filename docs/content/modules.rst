Python API
===========

Import sincei as:

.. code-block:: python

   import sincei

The following modules are available for use directly in python.

.. toctree::
   :maxdepth: 2
   :hidden:

   modules/ExponentialFamily
   modules/GLMPCA
   modules/TopicModels
   modules/multimodalClustering
   modules/ReadCounter
   modules/RegionQuery
   modules/WriteBedGraph
   modules/Utilities


* :doc:`modules/ReadCounter`
   Compute the coverage from bam files in genomic bins or regions (specified in BED or GTF files) and
   return a `cell x feature` count matrix. Can used for genome-wide coverage or for specified regions.

* :doc:`modules/ExponentialFamily`
   Implementation of many exponential family distributions to be used with the :doc:`modules/GLMPCA`.
   It can be easily extended to include other distributions.

* :doc:`modules/GLMPCA`
   Compute the generalized-linear model principal components (GLM-PCs) of a dataset, given an
   exponential distribution chosen based on prior knowledge.

* :doc:`modules/TopicModels`
   Compute LSA (Latent Semantic Analysis) or LDA (Latent Dirichlet Allocation) for an `AnnData`
   object and return the `cell x topic` matrix.

* :doc:`modules/multimodalClustering`
   Performs multi-graph clustering for cells with multi-modal data stored in a `MuData` object and
   return clusters and UMAP embeddings. 

* :doc:`modules/RegionQuery`
   Get overlaps between the features in an AnnData and regions in a BED or GTF file.

* :doc:`modules/WriteBedGraph`
   Write a bedgraph or bigwig file from bam files. Can be used for genome-wide coverage or for specified
   regions.

* :doc:`modules/Utilities`
   Utility functions used across the package.
