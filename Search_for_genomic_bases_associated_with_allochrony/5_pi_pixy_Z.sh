#!/bin/bash


gvcf="gVCF.portugal.vcf.gz"

bcftools -S SP_WP_males.list -r chrZ ${gvcf} -Oz > Z.SPWP.males.gVCF.vcf.gz
tabix -p vcf Z.SPWP.males.gVCF.vcf.gz

# Exécuter la commande pixy
pixy --stats pi \
     --vcf Z.SPWP.males.gVCF.vcf.gz \
     --populations SP_WP_males_pi.list \
     --window_size 2000 \
     --n_cores 6 \
     --output_folder calcul_pi_Z_SP_WP_males
