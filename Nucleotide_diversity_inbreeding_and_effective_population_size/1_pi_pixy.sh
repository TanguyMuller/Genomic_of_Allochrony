#!/bin/bash
# Script to calculate Pi for each chromosome

# gCVF file
gvcf="gVCF.portugal.vcf.gz"

# Extract the chromosome name
chr=$(basename "$gvcf" | cut -d'.' -f2)

# Create a folder for each chromosome
mkdir -p "$chr"

# Index the VCF file
tabix -f -p vcf "$gvcf"

# Run pixy to calculate nucleotide diversity (Pi)
# Note: population.portugal.list must be in 'ind<TAB>pop' format
pixy --stats pi \
   --vcf "$gvcf" \
   --populations population.portugal.list \
   --window_size 100000 \
   --n_cores 6 \
   --output_folder "$chr"
