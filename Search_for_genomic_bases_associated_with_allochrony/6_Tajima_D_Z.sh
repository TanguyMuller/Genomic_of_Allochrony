#!/bin/bash
# Script to compute Tajima's D on chromosome Z for LSP and LWP males
# Applies a 5% MAF filter and uses vcftools
# Usage: ./compute_TajimaD_chrZ.sh <window_size>

size=$1

# Extract samples for each population
bcftools view -S LSP_males.list Z.ind.portugal.vcf.gz -Oz -o Z_LSP.vcf.gz
bcftools view -S LWP_males.list Z.ind.portugal.vcf.gz -Oz -o Z_LWP.vcf.gz

# Compute Tajima's D with MAF ≥ 0.05
vcftools --gzvcf Z_LSP.vcf.gz --maf 0.05 --TajimaD $size --out chrZ_LSP_MAF005_${size}
vcftools --gzvcf Z_LWP.vcf.gz --maf 0.05 --TajimaD $size --out chrZ_LWP_MAF005_${size}
