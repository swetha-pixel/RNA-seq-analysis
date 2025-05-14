# RNA-seq-analysis
This RNA-seq analysis pipeline takes you from raw sequencing data to differential gene expression results. Follow these steps to process your data and identify differentially expressed genes:

Pre-Alignment QC: Perform quality control on raw sequencing data using FastQC.

Trimming : Trim adapters and low-quality bases from the reads if needed.

Alignment to Reference Genome: Align the reads to a reference genome using STAR.

Post-Alignment QC: Sort and index the BAM files, and check alignment quality using Samtools.

Transcript Quantification: Count the number of reads mapped to each gene using featureCounts.

Differential Gene Expression (DGE) Analysis: Perform DGE analysis using DESeq2 or similar packages in R to identify genes that are differentially expressed.

Visualization: Visualize results with plots such as volcano plots, heatmaps, and PCA to explore gene expression patterns.



This pipeline uses Mamba as the environment manager, though Conda can also be used.

## 1. Environment Setup (Mamba)

Create and activate a new environment:

```bash
mamba create -n rna_env python=3.8
mamba activate rna_env
```
Install required tools:

```bash
mamba install -c bioconda sra-tools fastqc trimmomatic star samtools subread
mamba install -c conda-forge r-base
```
## 2.Download the RNA-seq data

1.Extract the SRA toolkit:
```bash
tar -xvzf sratoolkit.3.0.0-ubuntu64.tar.gz
```
2.Download SRA data

```bash
fasterq-dump SRRXXXXX -S -p
```

## 3. Quality Control

Workflow of Quality control  looks like this:

1. Download FASTQ (fasterq-dump)
        ↓
2. Run FastQC on raw data  ← 🧪 Quality Check (pre-trim)
        ↓
3. Trim with trimmomatic  ← ✂️ Cleaning
        ↓
4. Run FastQC again on trimmed data ← 🧪 Quality Check (post-trim)
        ↓
5. Align / quantify / downstream analysis


▶️ Run FastQC

```bash
fastqc Sample_R1.fastq Sample_R2.fastq
```
✂️Trim with Trimmomatic

```bash
trimmomatic PE -threads 4 \
  Sample_R1.fastq Sample_R2.fastq \
  Sample_R1_paired.fastq Sample_R1_unpaired.fastq \
   Sample_R1_paired.fastq Sample_R2_unpaired.fastq \
  ILLUMINACLIP:/home/swetha/miniforge3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
```
## 4. Genome Indexing (STAR)

This step generates the genome index required for aligning RNA-seq reads using STAR. It uses the reference genome FASTA and GTF annotation file.

```bash
#!/bin/bash
# Activate your conda environment
mamba activate staril
# Create directory to store the STAR genome index
mkdir -p GenomeIndex
# Run STAR genome indexing
STAR \
  --runThreadN 4 \  # Adjust based on your CPU
  --runMode genomeGenerate \
  --genomeDir GenomeIndex/ \
  --genomeFastaFiles ./bombyx_mori.fa \
  --sjdbGTFfile ./bombyx_mori.gtf \
  --sjdbOverhang 99
```

## 5.RNA-seq Read Alignment (STAR)

Aligns quality-filtered RNA-seq reads to the indexed reference genome.
Generates sorted BAM files (by coordinate)

```bash
#!/bin/bash

# This script loops through multiple samples and performs alignment using STAR

# List of sample names (SRR IDs for each paired-end dataset)
samples=("SRRXXXXXXX4" "SRRXXXXXXX1")

# Define the temporary directory (relative or customizable path)
tmp_dir="./tmp/"

# Input and output base directories (relative or configurable)
input_dir="./data/trimmed/"
output_dir="./results/alignment/"
genome_dir="./GenomeIndex/"

# Create output and tmp directories if not present
mkdir -p "$output_dir" "$tmp_dir"

# Loop through each sample
for sample in "${samples[@]}"; do
  echo "Processing sample: $sample"

  # Define input files
  r1_file="${input_dir}${sample}_1.fastq.gz"
  r2_file="${input_dir}${sample}_2.fastq.gz"

  # Check if input files exist
  if [[ ! -f "$r1_file" || ! -f "$r2_file" ]]; then
    echo "Error: Input files for sample $sample not found."
    continue
  fi

  # Run STAR for alignment
  STAR --runMode alignReads \
    --runThreadN 4/6 \
    --readFilesIn "$r1_file" "$r2_file" \ #input SRR sample files
    --readFilesCommand zcat \
    --genomeDir "$genome_dir" \ #genome file path
    --outFileNamePrefix "${output_dir}${sample}_" \
    --outSAMtype BAM SortedByCoordinate \ #type of output
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --outTmpDir "${tmp_dir}${sample}_tmp/"

  # Check if STAR executed successfully
  if [[ $? -ne 0 ]]; then
    echo "Error: STAR alignment failed for sample $sample."
  else
    echo "Successfully processed sample: $sample"
  fi

done
```
## 6.Sorting BAM Files by Read Name

Re-sorts the alignment BAM files by read name using `samtools`. This is useful for tools that need the reads in read name order, such as certain RNA-seq quantifiers or visualization tools.

- Input: `_Aligned.sortedByCoord.out.bam` from STAR
- Output: `_sorted.bam` (sorted by read name)

```bash
#!/bin/bash

# This script sorts STAR-aligned BAM files by read name using samtools

# List of sample IDs
samples=("SRRXXXXXXX4" "SRRXXXXXXX1")

# Input and output directory paths (customizable or relative)
input_dir="./results/alignment/"
output_dir="./results/sorted_bam/"

# Check if the output directory exists; if not, create it
if [[ ! -d "$output_dir" ]]; then
  mkdir -p "$output_dir"
  echo "Created output directory: $output_dir"
fi

# Loop through each sample
for sample in "${samples[@]}"; do
  echo "Processing sample: $sample"

  input_bam="${input_dir}${sample}_Aligned.sortedByCoord.out.bam"
  output_bam="${output_dir}${sample}_sorted.bam"

  # Check if input BAM file exists
  if [[ ! -f "$input_bam" ]]; then
    echo "Error: Input BAM file for $sample not found. Skipping."
    continue
  fi

  # Use samtools to sort BAM file by read name
  samtools sort \
    -@ 4 \  # Number of threads (can be changed to 6)
    -n \
    -O BAM \
    "$input_bam" \# path of the input files
    -o "$output_bam"

  # Check if samtools executed successfully
  if [[ $? -ne 0 ]]; then
    echo "Error: samtools sorting failed for sample $sample."
  else
    echo "Successfully sorted BAM for sample: $sample."
  fi
done
```
## 7.Gene Quantification using featureCounts

This step uses `featureCounts` to count the number of reads mapping to genes, based on the aligned and sorted BAM files. The counts are used for downstream differential expression analysis.

- **Input:** Sorted BAM files (generated in the previous step)
- **GTF file:** Gene annotation file in GTF format
- **Output:** 
  - `featurecounts.txt` (Gene count table)
  - `featurecounts.txt.summary` (Summary of the quantification process)

```bash
#!/bin/bash
# genes_quantification_script.sh

# Activate Conda Environment (if required)
# source activate featurecounts

# Create the directory for storing counts if it doesn't exist
mkdir -p ./results/counts/

# Generate Counts Using FeatureCounts for paired-end reads
# Paired-end BAM files, sorted by name
featureCounts \
  -T 4 \                                      # Use 4 threads
  -p \                                       # Paired-end mode
  -S 2 \                                     # Strand-specific (reversely stranded)
  -a ./GenomeIndex/annotation.gtf \          # Path to GTF annotation
  -o ./results/counts/FeatureCounts.txt \    # Output file
  ./results/sorted_bam/*_sorted.bam          # Input BAM files
```
## 8. Differential Expression Analysis in R

All the following code should be executed in R.

This section outlines the steps for differential gene expression (DGE) analysis using DESeq2, starting from count data to quality control and functional interpretation.

📌 Overview of Steps:
Import Count Data

Load Data into DESeq2

Perform Differential Gene Expression Analysis

Generate QC Plots

(Optional) Functional Analysis – Enrichment of biological pathways

📊 Input Data
The primary input for DESeq2 is the raw count matrix (FeatureCounts_Mod.txt) with:

Row names as gene IDs
Columns representing sample-wise raw counts

# Install BiocManager if not already installed
install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install("SummarizedExperiment")
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
BiocManager::install("apeglm")
BiocManager::install("genefilter")
BiocManager::install("ggplot2", force = TRUE)

# Install CRAN packages
install.packages("gtable")
install.packages("stringi")
install.packages("tidyverse")
install.packages("pheatmap")
install.packages("PoiClaClu")
install.packages("readxl")



```bash
# Set working directory

setwd("C:/Users/Nagaraju/OneDrive/Desktop/Dokumen")
getwd()

#Load metadata

metadata <- read_excel("C:/Users/Nagaraju/OneDrive/Desktop/Dokumen/metadata.xlsx")
metadata <- as.data.frame(metadata)

#Clean metadata

metadata$Sample_ID <- trimws(metadata$Sample_ID)
metadata$Sample_ID <- toupper(metadata$Sample_ID)
metadata$Sample_ID <- as.factor(metadata$Sample_ID)
metadata$Condition <- as.factor(metadata$Condition)

#Check factor levels

levels(metadata$Sample_ID)
levels(metadata$Condition)

#Load raw gene expression data

counts <- read_excel("read_count.xlsx")
counts <- as.data.frame(counts)
rownames(counts) <- counts[[1]]
counts <- counts[,-1]
counts <- as.matrix(counts)

#Clean column names

colnames(counts) <- trimws(colnames(counts))
colnames(counts) <- toupper(colnames(counts))

#Create DESeq2 object

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ Condition)

#Normalize the data

dds <- DESeq(dds)

#Variance Stabilizing Transformation (for PCA & heatmap)

rld <- vst(dds)
head(assay(rld))

# Sample-to-sample distances

sampledists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampledists)

# Add condition labels to heatmap

rownames(sampleDistMatrix) <- paste(rld$Condition)
colnames(sampleDistMatrix) <- NULL

#Draw heatmap

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampledists,
         clustering_distance_cols = sampledists,
         col = colors)

#PCA Plot

plotPCA(rld, intgroup = c("Condition"))

#Ensure proper reference level

dds$Condition <- as.factor(dds$Condition)
dds$Condition <- relevel(dds$Condition, ref = "Control")

#Differential expression analysis

dds <- DESeq(dds)
res <- results(dds)

#Check metadata about results

mcols(res, use.names = TRUE)
summary(res)

# Filter significant genes

ressig <- subset(res, pvalue < 0.05)

# Inspect top genes

head(ressig[order(ressig$log2FoldChange), ])
head(ressig[order(-ressig$log2FoldChange), ])

# Plot expression of top gene

topGene <- rownames(res)[which.max(res$log2FoldChange)]
plotCounts(dds, gene = topGene, intgroup = c("Condition"))

#Optional: Get exact min and max programmatically(xlim and ylim)

min_log2FC <- min(ressig$log2FoldChange, na.rm = TRUE)
max_log2FC <- max(ressig$log2FoldChange, na.rm = TRUE)
cat("Range of log2 Fold Change:", min_log2FC, "to", max_log2FC, "\n")


#Volcano plot

plot(EnhancedVolcano(res,
                     lab = rownames(res),
                     x = 'log2FoldChange',
                     y = 'pvalue',
                     pCutoff = 1e-2,
                     FCcutoff = 1.00,
                     xlim = c(-27, 27),
                     ylim = c(0, -log10(min(res$pvalue, na.rm=TRUE))),
                     pointSize = 1.3,
                     labSize = 2.6,
                     max.overlaps = 100,                    # increase label capacity
                     title = "Control vs Infected",
                     subtitle = "Differential expression analysis",
                     caption = "Log2FC cutoff = 1; p-value cutoff = 1e-2",
                     legendPosition = "right",
                     colAlpha = 0.6,
                     drawConnectors = TRUE,
                     hline = -log10(1e-2),                 # horizontal line at pCutoff
                     widthConnectors = 0.5))


#Save results to file

resord <- as.data.frame(res)
finaltable <- resord[order(resord$pvalue), ]
write.csv(finaltable, file = "finaltable.csv", row.names = TRUE, quote = TRUE)
```
