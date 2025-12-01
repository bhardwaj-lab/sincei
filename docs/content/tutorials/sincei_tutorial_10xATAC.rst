Analysis of 10x genomics scATAC-seq data using sincei
=====================================================

The data used here is produced in our :doc:`10x genomics multiome tutorial <sincei_tutorial_10x>`.

We will need a blacklist file to filter out reads aligned to known problematic regions of the human
genome. This file and blacklists for other genome assemblies can be downloaded from the `Boyle lab
ENCODE blacklist repository <https://github.com/Boyle-Lab/Blacklist>`__.

.. code:: bash

   wget -O hg38-blacklist.v2.bed.gz https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz

We will use the BAM files, barcodes and peaks.bed file obtained from the **cellranger-arc**
workflow. Alternatively, the BAM file (tagged with cell barcodes) and peaks.bed can be obtained from
a custom mapping/peak calling workflow, and the list of filtered cell barcodes can be obtained using
sincei (see the :doc:`parent tutorial <sincei_tutorial_10x>` for explanation).

Define common bash variables:

.. code:: bash

   # create dir
   mkdir sincei_output/atac

   # save as bash variables
   blacklist=hg38-blacklist.v2.bed.gz

2. Quality control - I (read-level)
-----------------------------------

In order to identify high quality cells for our analysis, we can use the read-level quality
statistics from :ref:`scFilterStats`. Low quality cells in this data can be identified using several
criteria, such as:

- high number of PCR duplicates (filtered using `--duplicateFilter`)
- high fraction of reads aligned to blacklisted regions (filtered using `--blacklist`)
- high fraction of reads with poor mapping quality (filtered using `--minMappingQuality`)
- vey high/low GC content of the aligned reads, indicating the reads were mostly aligned to
  low-complexity regions (filtered using `--GCcontentFilter`)
- high level of secondary/supplementary alignments (filtered using `--samFlagExclude/Include`)

.. code:: bash


   for rep in rep1 rep2
   do
      dir=cellranger_output_${rep}/outs/
      bamfile=${dir}/atac_possorted_bam.bam
      barcodes=${dir}/filtered_feature_bc_matrix/barcodes.tsv.gz # from sincei or cellranger output
      zcat ${barcodes} > ${dir}/filtered_barcodes.txt

      scFilterStats -p 20 \
         --region chr1 \
         --GCcontentFilter '0.2,0.8' \
         --minMappingQuality 10 \
         --samFlagInclude 64 \
         --samFlagExclude 256 \
         --samFlagExclude 2048 \
         --barcodes ${dir}/filtered_barcodes.txt \
         --cellTag CB \
         -bl ${blacklist}  \
         --label atac_${rep} \
         -o sincei_output/atac/scFilterStats_output_${rep}.txt \
         -b ${bamfile}
   done

:ref:`scFilterStats` summarizes these outputs as a table, which can then be visualized using the
`MultiQC tool <https://multiqc.info/docs/>`__, to select appropriate list of cells to include for
counting.

.. code:: bash

   multiqc sincei_output/atac/ # results in multiqc_report.html

3. Signal aggregation (counting reads)
--------------------------------------

Below we aggregate the signal from single-cells, using the ``atac_peaks.bed`` file, which contains
regions of local signal enrichment detected using the ``cellranger-arc`` workflow.

Since we get a seperate peaks file from each replicate, we can run sincei on the union of the peaks
detected on each sample.

If needed, we can use the same parameters as in :ref:`scFilterStats` to count only high quality reads
from our whitelist of barcodes. We avoid counting reads in blacklisted regions of the human genome
(--blacklist).

.. code:: bash

   ## merge intervals from 2 peaks bed files
   for f in cellranger_output_rep*/outs/atac_peaks.bed
   do
      awk '{if(NR>51) {print $0}}' $f >> repmerged.bed
   done

   sort -k1,1 -k2n,2 repmerged.bed | bedtools merge -i - > atac_peaks_merged.bed

   # count reads on merged bed file
   for rep in rep1 rep2
   do
      dir=cellranger_output_${rep}/outs/
      bamfile=${dir}/atac_possorted_bam.bam
      barcodes=${dir}/filtered_barcodes.txt

      scCountReads features -p 20 \
         --BED atac_peaks_merged.bed --cellTag CB \
         --region chr1 \
         --minMappingQuality 10 \
         --samFlagInclude 64 \
         --samFlagExclude 2048 \
         --extendReads \
         -bl ${blacklist} \
         -bc ${barcodes} \
         --cellTag CB \
         -o sincei_output/atac/scCounts_atac_peaks_${rep} \
         --label atac_${rep} \
         -b ${bamfile}
   done
   # Number of bins found: 18527

4. Quality control - II (count-level)
-------------------------------------

After counting, it is recommended to perform QC on these counts, in order to filter regions and
cells that have low counts, or have low enrichment of counts. Even though we already performed
read-level QC before, the counts distribution on our specified regions (bins/genes/peaks) could be
different from the whole-genome stats.

We can run :ref:`scCountQC` on the count data to get various statistics at region and cell level.
Running this tool with the ``--describe`` flag lists the metrics that can be used to filter
cells/regions.

.. code:: bash

   # list the metrics we can use to filter cells/regions
   for rep in rep1 rep2
   do
      scCountQC -i sincei_output/atac/scCounts_atac_peaks_${rep}.h5ad --describe
   done

The tool :ref:`scCountQC` can be used for count-level QC and filtering of count data. With the
``--outMetrics`` option, the tool outputs the count statistics at region and cell level (labelled as
``<prefix>.regions.tsv`` and ``<prefix>.cells.tsv``). Just like :ref:`scFilterStats`, these outputs
can then be visualized using the `MultiQC tool <https://multiqc.info/docs/>`__, to select
appropriate metrics to filter out the unwanted cells/regions.

.. code:: bash

   # export the single-cell level matrices
   for rep in rep1 rep2
   do
      scCountQC -i sincei_output/atac/scCounts_atac_peaks_${rep}.h5ad -om sincei_output/atac/countqc_atac_peaks_${rep}
   done

   # visualize output using multiQC
   multiqc sincei_output/atac/ # see results in multiqc_report.html

In this example, we detect ~18.5k peaks in ~13.5K cells.

Below, we perform some basic filtering using :ref:`scCountQC`. We exclude the cells with <500 and
>10000 detected peaks (``--filterRegionArgs``). We also exclude peaks detected in too few cells
(<100) or in >90% of cells (``--filterCellArgs``).

.. code:: bash

   for rep in rep1 rep2
   do
      scCountQC -i sincei_output/atac/scCounts_atac_peaks_${rep}.h5ad \
         -o sincei_output/atac/scCounts_atac_peaks_filtered_${rep} \
         --filterRegionArgs "n_cells_by_counts: 100, 5500" \
         --filterCellArgs "n_genes_by_counts: 500, 10000"
   done

   ## rep 1
   #Applying filters
   # Cells post-filtering: 6153
   # Features post-filtering: 13502
   ## rep 2
   #Applying filters
   # Cells post-filtering: 6109
   # Features post-filtering: 13911

5. Combine counts for the 2 replicates
--------------------------------------

While it is recommended to perform count QC separately for each replicate, the counts can now be
combined into one file for more convenient for downstream analysis. We provide a tool
:ref:`scCombineCounts`, which can concatenate counts for cells based on common features (in
``multi-sample`` mode). It can also be used to concatenate different modalities based on a common
set of cells (in ``multi-modal`` mode).

Concatenating the filtered cells for the 2 replicates results in a total of ~12K cells.

.. code:: bash

   scCombineCounts \
      -i sincei_output/atac/scCounts_atac_peaks_filtered_rep1.h5ad \
      sincei_output/atac/scCounts_atac_peaks_filtered_rep2.h5ad \
      -o sincei_output/atac/scCounts_atac_peaks_filtered.merged.h5ad \
      --method multi-sample \
      --labels rep1 rep2
   # Combined cells: 12262
   # Combined features: 13249

5. Dimensionality reduction and clustering
------------------------------------------

The tool :ref:`scClusterCells` provides a range of options to reduce the dimensionality of our count
data, while preserving biological signal. This can be specified with ``--method`` option. Below, we
will use a topic modeling method called Latent Semantic Analysis (LSA) to reduce the dimensionality
of our data to 30 topics (the default). The tool then performs Leiden clustering, and presents a
UMAP (dimensionality reduction to 2 dimensions) plot of the output (``--outFileUMAP`` option). This
option also creates a tsv file with the UMAP coordinates and assigned cluster for each cell in our
data.

.. code:: bash

   scClusterCells -i sincei_output/atac/scCounts_atac_peaks_filtered.merged.h5ad \
      -m LSA \
      -n 30 \
      --clusterResolution 1 \
      -op sincei_output/atac/scClusterCells_UMAP.png \
      -o sincei_output/atac/scCounts_atac_peaks_clustered.h5ad
   # Coherence Score:  -1.83
   # also produces the tsv file "sincei_output/scClusterCells_UMAP.tsv"

(optional) Confirmation of clustering using metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below, we load this data in R and compare it to the cell metadata provided with our files to verify
that our clustering separates celltypes in a biologically meaningful way.

We can color our UMAP output from :ref:`scClusterCells` with the cell-type information from `Persad et.
al. (2023) <https://www.nature.com/articles/s41587-023-01716-9>`__, that we provide `here
<https://figshare.com/articles/dataset/10x_multiome_test_data_package/29424470>`__.

.. code:: r

   library(dplyr)
   library(magrittr)
   library(ggplot2)
   library(patchwork)

   umap <- read.delim("sincei_output/atac/scClusterCells_UMAP.tsv")
   meta <- read.csv("metadata_cd34_atac.csv", row.names = 1)
   umap$celltype <- meta[gsub("rep1_", "", umap$Cell_ID), "celltype"]

   # keep only cells with published labels
   umap %<>% filter(!is.na(celltype))

   ## make plots
   df_center <- group_by(umap, cluster) %>% summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))
   df_center2 <- group_by(umap, celltype) %>% summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

   # colors for metadata (12 celltypes)
   col_pallete <- c("#808080", RColorBrewer::brewer.pal(8, "Paired"))
   names(col_pallete) <- unique(umap$celltype) # grey is for NA

   # colors for sincei UMAP (9 clusters)
   colors_cluster <- c("#315e66", "#808080", "#315e66", RColorBrewer::brewer.pal(12, "Paired"))
   names(colors_cluster) <- unique(umap$cluster)


   p1 <- umap %>%
     ggplot(., aes(UMAP1, UMAP2, color=factor(cluster), label=cluster)) +
     geom_point() + geom_label(data = df_center, aes(UMAP1, UMAP2)) +
     scale_color_manual(values = colors_cluster) + theme_minimal(base_size = 12) +
     theme(legend.position = "none") + #labs(x="UMAP1", y="UMAP2") +
     ggtitle("sincei clusters (LSA + Leiden)")

   p2 <- umap %>% filter(!is.na(celltype)) %>%
     ggplot(., aes(UMAP1, UMAP2, color=factor(celltype), label=celltype)) +
       geom_point() + geom_label(data = df_center2, aes(UMAP1, UMAP2)) +
       scale_color_manual(values = col_pallete) +
       labs(color="Cluster") + theme_minimal(base_size = 12) +
       theme(legend.position = "none") +
     ggtitle("Published Cell Types")

   pl <- p1 + p2
   pl
   ggsave(plot=pl, "sincei_output/atac/UMAP_compared_withOrig.png", dpi=300, width = 11, height = 6)


.. image:: ./../images/UMAP_compared_withOrig_10xATAC.png
    :height: 800px
    :width: 1600 px
    :scale: 50 %

The figure above shows that we can easily replicate the expected cell-type results from the scATAC
data using **sincei**. This was done using only 1/23th of original data (chromosome 1) and basic
pre-processing steps, therefore the results should only improve with full data, better cell/region
filtering and optimizing the analysis parameters.

6. Creating bigwigs and visualizing signal on IGV
-------------------------------------------------

For further exploration of data, it can be useful to create pseudo-bulk coverage files (bigwigs)
that aggregate the signal across cells in our each of our clusters. The tool :ref:`scBulkCoverage`
takes the clustering information `.tsv` file produced by :ref:`scClusterCells`, along with the
corresponding BAM files, and aggregates the signal to create these bigwigs.

The parameters here are same as other sincei tools that work on BAM files, except that we can output
normalized bulk signal (specified using `--normalizeUsing` option) . Below, we produce
CPM-normalized bigwigs at 1kb resolution.

.. code:: bash

   scBulkCoverage -p 20 \
   --cellTag CB \
   --normalizeUsing CPM \
   --binSize 1000 \
   --minMappingQuality 10 \
   --samFlagInclude 64 \
   --samFlagExclude 2048 \
   --duplicateFilter 'start_bc_umi' \
   --extendReads \
   -b cellranger_output_rep1/outs/atac_possorted_bam.bam \
   cellranger_output_rep2/outs/atac_possorted_bam.bam \
   --labels rep1_atac_rep1 rep2_atac_rep2 \
   -i sincei_output/atac/scClusterCells_UMAP.tsv \
   -o sincei_output/atac/sincei_cluster
   # creates 9 files with names "sincei_cluster_<X>.bw" where X is 0, 1... 9

We can now inspect our bigwigs on `IGV <https://software.broadinstitute.org/software/igv/>`__. We
can clearly see some regions with cell-type specific signal, such as the markers described in the
original manucscript: TAL1 (erythroid), MPO (myeloid) and IRF (dendritic) marker genes.

.. image:: ./../images/igv_snapshot_10xATAC.png
   :height: 500px
   :width: 6000 px
   :scale: 50 %
