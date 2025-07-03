#!/bin/bash
# Script to calculate linkage disequilibrium (LD) decay using LSP individuals
# Input: LSP-specific VCF files (subset from ind.portugal.freebayes_filt_maf005.vcf.gz for autosomes 
#        and ind.portugal.freebayes.chrZ_filt_maf005.vcf.gz for chromosome Z)

# Input arguments
VCF_FILE=$1
CHR=$2

# Set the chromosome parameter
if [[ $CHR == "Auto" ]]; then
    CHR_SET='--chr-set 49'  # 49 autosomal chromosomes
    CHR_TARGET='chr1'
elif [[ $CHR == "chrZ" ]]; then
    CHR_SET='--chr-set 1'   # 1 Z chromosome
    CHR_TARGET='chrZ'
else
    echo "Error: Invalid chromosome parameter. Use 'Auto' or 'chrZ'."
    exit 1
fi

# Calculate LD decay on LSP individuals 
./plink --vcf $VCF_FILE --double-id --allow-extra-chr \
$CHR_SET \
--chr $CHR_TARGET \
--thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 \
--ld-window-r2 0 \
--make-bed --out ld_decay/ld_decay_thin0.1_${CHR}

# Run Python script for identifying the LD decay in LSP
./ld_decay_calc.py -i ld_decay/ld_decay_thin0.1_${CHR}.ld.gz -o ld_decay/ld_decay_thin0.1_${CHR}

echo "End LD decay analysis"

### Commands to run 2_ld_decay_LSP.sh
#bash run 2_ld_decay_LSP.freebayes_filt_maf005.sh LSP.vcf.gz Auto
#bash run 2_ld_decay_LSP.freebayes_filt_maf005.sh LSP.vcf.gz chrZ
