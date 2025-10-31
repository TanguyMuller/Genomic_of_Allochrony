#! /bin/bash
# Script for create ped and map file for GONe analyses
# and for currentNe2

mkdir -p vcf
pop=$1

VCF="${pop}_ind_masked_HWE.vcf.gz"
plink --vcf ${VCF} --recode --out vcf/${pop}_prune_mono --chr-set 49 --allow-extra-chr --double-id --mac 1
awk '{print $1,NR,$3,$4}' vcf/${pop}.pruned_mono.map > vcf/tmp.${pop}.pruned_mono.MAF.map
mv vcf/${pop}.pruned_mono.map vcf/${pop}.pruned_mono.map.unchanged
mv vcf/tmp.${pop}.pruned_mono.MAF.map vcf/${pop}.pruned_mono.map
