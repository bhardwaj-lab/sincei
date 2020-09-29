# sincei

Single-Cell Informatics: A user-friendly toolkit for QC, counting, clustering and plotting of single-cell (epi)genomics data.


### Currently supported functions

1. scFilterStats: Evaluate cells from a BAM file based on various QC metrics. Output: tsv file

2. scCountReads: Count reads/UMIs on genes, bins, exons etc, from a BAM file. Output: mtx file / hdf5 object (for other tools).

3. scClusterCells: Cluster cells using LSA algorithm and make UMAP. Input: hdf5 object (from scCountReads), Output: tsv file (with UMAP + clusters), UMAP plot.

4. scBulkCoverage: Create pseudo-bulk bigWig files from BAM files based on cluster assignments (from scClusterCells). Input: BAM file + tsv file (from scClusterCells), Output: bigWig files per cluster

### Other tool ideas

5. scNucleoFilter: Evaluate nucleosome banding patterns per cell (for scATAC/scChIC/scCutNRun protocols). Input: BAM-files, Output: tsv-file (with bandwidth values per cell) + plot (to detect low quality/overdigested cells).

6. scFingerPrint: Evaluate cells based on read enrichment distribution in genomic bins. Identify clusters of cells with high/low enriched genomic fraction based on histone mark of interest (for scChIC/scCutNRun protocols). Input: BAM-files, Output: tsv-files + UMAP plot.

7. scPlotEnrichment: Plot the fraction of read enrichment per cell on features from a BED/GTF file. Input: BAM-file, Output: Enrichment plot + tsv file.

8. scPlotRegion: Make a plot of read counts per cell in ordered genomic windows. Input: hdf5 file (from scCountReads), Output: plot of counts (cell\*regions)

9. scFindMarkers: Find most variable regions/genes per cell cluster based on a statistical test. Input: hdf5 file, tsv file (with clusters), Output: tsv file, with marker genes per cluster + pvalues
