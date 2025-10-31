#! /bin/bash
# Script for create ped and map file for GONe analyses

pop=$1

VCF="${pop}_ind_masked_HO.cvf.gz"
#bcftools view -T ^list.position.N.txt vcf/${VCF} -Oz > vcf/masked.${VCF}
plink --vcf vcf/masked.${VCF} --recode --out vcf/masked.chrA.${pop}.pruned_mono.MAF --chr-set 49 --allow-extra-chr --double-id --mac 1
awk '{print $1,NR,$3,$4}' vcf/masked.chrA.${pop}.pruned_mono.MAF.map > vcf/tmp.masked.chrA.${pop}.pruned_mono.MAF.map
mv vcf/masked.chrA.${pop}.pruned_mono.MAF.map vcf/masked.chrA.${pop}.pruned_mono.MAF.map.unchanged
mv vcf/tmp.masked.chrA.${pop}.pruned_mono.MAF.map vcf/masked.chrA.${pop}.pruned_mono.MAF.map 
