#!/bin/bash 
# Script to filter SNPs without 10 reads per pool
# and estimate the allele frequency from pools
# for RZooRoH analyses.
#
# Command
# bash 4_freq_alleliq_pool.sh ind.vcf.gz Auto

vcf=$1
path_list="/your/path/to/list"
path_script="/your/path/to/script"


# Set directory for RZooRoH
mkdir -p freq_alleliq
mkdir -p results
mkdir -p GL

# Last filt on min_cov_per_pool=10 
# Create freq alleliq on pool for each population
Rscript ${path_script}/freq_alleliq_pool_exp.R ${dataset}

# Keep only autosomes SNPs
grep -v chrZ ${path_list}/positions_${base_name}_10covperpool.txt > tmp.list
mv tmp.list ${path_list}/positions_${base_name}_10covperpool.txt

bcftools view -T ${path_list}/positions_${base_name}_10covperpool.txt -Oz ${path_vcf}/${base_name}.vcf.gz > ${path_vcf}/${base_name}_10covperpool.vcf.gz
bcftools index ${path_vcf}/${base_name}_10covperpool.vcf.gz
#rm ${path_vcf}/${dataset}*

#if [[ "$dataset" != "chr49_filt_0pct_missing_dp_filtered_pool.vcf.gz" ]]; then
#    rm -f ${path_vcf}/${base_name}.vcf.gz*
#fi

gunzip -c ${path_vcf}/${base_name}_10covperpool.vcf.gz | grep -c ^c
