#!/bin/bash 
# Script to filter SNPs with fewer than 10 reads per pool
# and to estimate allele frequencies from pool data
# for RZooRoH analyses.
#
# Command
# bash 4_freq_alleliq_pool.sh ind.vcf.gz 

vcf=$1
path_list="/your/path/to/list"
path_script="/your/path/to/script"
base_name=$(basename "$vcf" .vcf.gz)

# Prepare output directories
mkdir -p freq_alleliq
mkdir -p results
mkdir -p GL

# Run the R script that computes allele frequencies per pool
echo "Estimating allele frequencies for ${vcf}..."
Rscript ${path_script}/freq_alleliq_pool_exp.R ${vcf}

# Remove sites with fewer than 10 reads per pool
echo "Filtering SNPs with at least 10 reads per pool..."
bcftools view -T ${path_list}/positions_${vcf}_10covperpool.txt -Oz ${vcf} > ${base_name}_10covperpool.vcf.gz
bcftools index ${base_name}_10covperpool.vcf.gz

# Number of SNPs retained for RZooRoH analyses
num_snps=$(gunzip -c "${base_name}_10covperpool.vcf.gz" | grep -c ^chr)
echo "Number of SNPs retained for RZooRoH analyses: ${num_snps}"

echo "Filtered VCF created: ${base_name}_10covperpool.vcf.gz"
echo "Process completed successfully"
