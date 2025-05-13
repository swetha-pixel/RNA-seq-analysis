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
