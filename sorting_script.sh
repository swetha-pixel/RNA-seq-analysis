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
