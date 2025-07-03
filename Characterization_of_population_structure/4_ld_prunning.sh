#!/bin/bash
# Script for pruning VCF by linkage disequilibrium (LD)
# Creates prune.in and prune.out files using PLINK from MAF-filtered VCF files

# Input arguments
FILE=$1
CHR=$2
base_name=$(basename "$FILE" .vcf.gz)

# Set chromosome parameters
if [[ $CHR == "Auto" ]]; then
    CHR_SET='--chr-set 49'  # 49 autosomal chromosomes
elif [[ $CHR == "chrZ" ]]; then
    CHR_SET='--chr-set 1'   # 1 Z chromosome
else
    echo "Error: Invalid chromosome parameter. Use 'Auto' or 'chrZ'."
    exit 1
fi

# Execute PLINK commands for individual-based VCF files
# Using 20kb window size based on LD decay analysis from previous scripts
if [[ $base_name == "ind.portugal.freebayes_filt_maf005" || $base_name == "ind.portugal.freebayes.chrZ_filt_maf005" ]]; then
    ./plink --vcf $FILE --indep-pairwise 20kb 1 0.1 --allow-extra-chr $CHR_SET --out $base_name
    awk '{split($1,tmp,":");{print tmp[1]"\t"tmp[2]}}' ${base_name}.prune.in > ${base_name}_good_positions.list
fi

# Create the pruned dataset keeping only unlinked SNPs
bcftools view -T ${base_name}_good_positions.list $FILE -Oz > ${base_name}_ldprune.vcf.gz

echo "Number of unlinked SNPs: $(bcftools view -H ${base_name}_ldprune.vcf.gz | wc -l)" 

### Commands to run 4_ld_pruning.sh
#bash 4_ld_pruning.sh ind.portugal.freebayes_filt_maf005.vcf.gz Auto
#bash 4_ld_pruning.sh ind.portugal.freebayes.chrZ_filt_maf005.vcf.gz chrZ
#bash 4_ld_pruning.sh all.portugal.freebayes_filt_maf005.vcf.gz Auto
#bash 4_ld_pruning.sh all.portugal.freebayes.chrZ_filt_maf005.vcf.gz chrZ
#bash 4_ld_pruning.sh pool.portugal.freebayes_filt_maf005.vcf.gz Auto
#bash 4_ld_pruning.sh pool.portugal.freebayes.chrZ_filt_maf005.vcf.gz chrZ
