Analysis of 10x genomics multiome data using sincei
===================================================

Below, we will demonstrate how we can use sincei to explore the
scRNA-seq and scATAC-seq data as part of the `10x multiome
protocol <https://www.10xgenomics.com/products/single-cell-multiome-atac-plus-gene-expression>`__.
The 10x multiome kit allows joint profiling of single-cell ATAC-seq and
RNA-seq from single-cells. Here, we will analyse these two data sets
separately. We will will use the dataset published with `Persad et.
al. (2023) <https://www.nature.com/articles/s41587-023-01716-9>`__,
which profiles CD34+ cells from human bone marrow.

1. Download and process the dataset
-----------------------------------

The raw fastq files were downloaded from
`GEO <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200046>`__
and processed using the standard 10x genomics `cellranger-arc
workflow <https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/algorithms/overview>`__.

.. code:: bash

    # Download using sra-tools
    prefetch --force all --max-size 10t -O GSM6005302_CD34_Rep1_Multiome_RNA_Homo_sapiens_RNA-Seq SRR18590099 SRR18590100 && \
    cd GSM6005302_CD34_Rep1_Multiome_RNA_Homo_sapiens_RNA-Seq && vdb-validate SRR18590099 SRR18590100 && \
    fasterq-dump -m 10000MB -b 100MB -c 1000MB --temp /<path>/<to>/tempdir \
    -e 20 -p --split-files --include-technical -O . SRR18590099 SRR18590100

    # Run cellranger-arc
    for rep in rep1 rep2
    do cellranger-arc count --disable-ui --reference \
      /<path>/<to>/cellranger_indicies/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/ \
      --localcores 20 --localmem 200 --id output_${rep} --libraries cellranger_samplesheet_${rep}.csv
    done

Below is the structure of samplesheet that was used for cellranger:

.. code:: bash

    fastqs,sample,library_type
    /path/to/10x_multiome/fastq,CD34_Rep1_Multiome_ATAC,Chromatin Accessibility
    /path/to/10x_multiome/fastq,CD34_Rep1_Multiome_RNA,Gene Expression


Below is the structure of the output directory from the workflow:

.. code:: bash

   <output_di>/outs:
   ├── analysis
   ├── atac_cut_sites.bigwig
   ├── atac_fragments.tsv.gz
   ├── atac_fragments.tsv.gz.tbi
   ├── atac_peak_annotation.tsv
   ├── atac_peaks.bed
   ├── atac_possorted_bam.bam
   ├── atac_possorted_bam.bam.bai
   ├── cloupe.cloupe
   ├── filtered_feature_bc_matrix
   ├── filtered_feature_bc_matrix.h5
   ├── gex_molecule_info.h5
   ├── gex_possorted_bam.bam
   ├── gex_possorted_bam.bam.bai
   ├── per_barcode_metrics.csv
   ├── raw_feature_bc_matrix
   ├── raw_feature_bc_matrix.h5
   ├── summary.csv
   └── web_summary.html

We will use the ``gex_possorted_bam.bam`` for gene-expression data and
``atac_possorted_bam.bam`` for chromatin accessibility analysis using
sincei. These files can also be produced as part of the
``cellranger count`` workflow for scRNA-seq or scATAC-seq data alone.
For convenience, we provide a subset of this data (only a small chunk of chromosome 2)
`here <https://figshare.com/articles/dataset/10x_multiome_test_data_package/29424470>`__

.. code:: bash

   mkdir 10x_multiome && wget -O 10x_multiome/10x_multiome_testdata.tar.gz https://figshare.com/ndownloader/files/55726430
   tar -xvzf 10x_multiome/10x_multiome_testdata.tar.gz ## releases 4 (indexed) bam files and 2 metadata files.

(optional) pre-filtering of barcodes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most of the cell barcodes from the droplet-based protocols (like 10x
genomics) do not contain cells. Therefore they have very low counts.
These must be filtered away at the beginning of the analysis. Although
the cellranger pipeline already provides a list of filtered barcodes,
sincei also allows you to extract per barcode count distributions,
indicating which barcodes should be removed. This can be done using the
``scFilerBarcodes`` tool.

.. code:: bash

   barcodes=737K-arc-v1.txt # cellranger-arc barcodes in this case
   for r in 1 2
   do
       bamfile=cellranger_output_rep${r}/outs/atac_possorted_bam.bam
       scFilterBarcodes -p 20 -b ${bamfile} -w ${barcodes} \
       -o sincei_output/atac_barcodes_rep${r}.tsv \
       --minCount 100 --minMappingQuality 10 --cellTag CB \
       --rankPlot sincei_output/barcode_rankplot_rep${r}.png
   done

The above example uses a whitelist of possible ATAC barcodes from
``cellrange-arc`` workflow. `See
here <https://kb.10xgenomics.com/hc/en-us/articles/360049105612-Barcode-translation-in-Cell-Ranger-ARC>`__
for more details. Providing a whitelist is optional in general, but
recommended for 10x genomics data.

The output file contains a list of filtered barcodes that contain counts
in atleast ``-mc`` regions of the genome. Unlike other tools with
similar options, sincei splits the data in 100kb bins and reports
whether or not a barcode has signal in those bins. This way, barcodes
with high counts, but present in only one genomic bin can also be
filtered out. In most cases, the output is same as the usual approach of
filtering by total counts. ``-rp`` produces the familiar ``knee-plot``
of the barcode counts.

2. scATAC-seq analysis
----------------------

Please follow :doc:`this tutorial <sincei_tutorial_10xATAC>` for further analysis of scATAC-seq samples from the above data.

3. scRNA-seq analysis
---------------------

Please follow :doc:`this tutorial <sincei_tutorial_10xRNA>` for further analysis of scRNA-seq samples from the above data.

Notes
------------

Currently, sincei doesn't provide a method for **doublet estimation and removal**, which is an important step in the analysis of droplet-based data.
Instead, we use simpler filters of min and max number of detected features per cell, which, to some extent mitigates this issue. However, this could
lead to some differences in results compared to the published data in used here. Despite this difference, the major published cell types can be separated
with sincei for both ATAC and RNA fraction of the data, as shown in the 2 tutorials above.
