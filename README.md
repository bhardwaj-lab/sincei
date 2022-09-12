# sincei

[![DOI](https://zenodo.org/badge/271841139.svg)](https://zenodo.org/badge/latestdoi/271841139)

Bhardwaj V. (2022) Single-Cell Informatics: A user-friendly toolkit for QC, counting, clustering and plotting of single-cell (epi)genomics data.


### Currently supported functions

sincei is a suite of command-line tools developed for a user-friendly analysis of single-cell sequencing data.
Version: 0.1

**Usage**: get the tool list with `sincei --help`

Each tool begins with the prefix sc<tool_name>, such as:

 $ scBulkCoverage -b file1.bam -g groupinfo.txt -o coverage

[ Tools for a typical single-cell analysis workflow ] (WIP: work in progress/not available yet)

    scFilterBarcodes        Identify and filter cell barcodes from BAM file (for droplet-based single-cell seq)
    scFilterStats           Produce per-cell statistics after filtering reads by user-defined criteria.
    scCountReads            Counts reads for each barcode on genomic bins or user-defined features.
    scCountQC               Perform quality control and filter the output of scCountReads.
    scCombineCounts         [WIP] Concatenate/merge the counts from different samples/batches or modalities
    scClusterCells          Perform dimensionality reduction and clustering on the output of scCountReads.
    scFindMarkers           [WIP] Find marker genes per group, given the output of scCountReads and a user-defined group.
    scBulkCoverage          Get pseudo-bulk coverage per group using a user-supplied cell->group mapping (output of scClusterCells).
    scFeaturePlot           [WIP] Plot the counts for a given feature on a UMAP or on a (IGV-style) genomic-track.
