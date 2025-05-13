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
