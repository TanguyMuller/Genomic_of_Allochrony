#!/bin/bash
# Script to calculate Pi for the Z chromosome for LSP and LWP males

# gCVF file
gvcf="gVCF.portugal.vcf.gz"

bcftools -S SP_WP_males.list -r chrZ ${gvcf} -Oz > Z.SPWP.males.gVCF.vcf.gz
tabix -p vcf Z.SPWP.males.gVCF.vcf.gz

# Run pixy to calculate nucleotide diversity (Pi)
# Note: SP_WP_males_pi.list must be in 'ind<TAB>pop' format
pixy --stats pi \
     --vcf Z.SPWP.males.gVCF.vcf.gz \
     --populations SP_WP_males_pi.list \
     --window_size 2000 \
     --n_cores 6 \
     --output_folder calcul_pi_Z_SP_WP_males
