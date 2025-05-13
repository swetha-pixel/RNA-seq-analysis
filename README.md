# RNA-seq-analysis
This RNA-seq analysis pipeline takes you from raw sequencing data to differential gene expression results. Follow these steps to process your data and identify differentially expressed genes:

Pre-Alignment QC: Perform quality control on raw sequencing data using FastQC.

Trimming (Optional): Trim adapters and low-quality bases from the reads if needed.

Alignment to Reference Genome: Align the reads to a reference genome using STAR.

Post-Alignment QC: Sort and index the BAM files, and check alignment quality using Samtools.

Transcript Quantification: Count the number of reads mapped to each gene using featureCounts.

Differential Gene Expression (DGE) Analysis: Perform DGE analysis using DESeq2 or similar packages in R to identify genes that are differentially expressed.

Visualization: Visualize results with plots such as volcano plots, heatmaps, and PCA to explore gene expression patterns.

#RNA-Seq Data Download and Environment Setup
This section will guide you through setting up the environment and downloading RNA-seq data using SRR IDs.

1. 1. Set up your environment 
This pipeline uses Mamba as the environment manager, though Conda can also be used.

## 🔧 Environment Setup (Mamba)

Create and activate a new environment:

```bash
mamba create -n rna_env python=3.8
mamba activate rna_env

Install required tools:
mamba install -c bioconda sra-tools fastqc trimmomatic star samtools subread
mamba install -c conda-forge r-base
mamba install -c bioconda bioconductor-deseq2
mamba install -c conda-forge r-tidyverse

