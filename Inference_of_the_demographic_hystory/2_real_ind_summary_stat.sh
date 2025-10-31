#! /bin/bash
# Script for estimate individual summary statistics on real data

# VCF with individual and pool data of LSP, LWP and FU
INPUT_VCF="demo.vcf.gz"
LIST="/path/to/your/list"

# Heterozygosity
vcftools --gzvcf $INPUT_VCF --hardy --keep ${LIST}/SP.list --out SP
awk 'NR>1' SP.hwe | awk '{split($3,el,"/"); if ((el[1]+el[2]+el[3])>0) print (el[2]/(el[1]+el[2]+el[3]))}' > SP.het.txt
vcftools --gzvcf $INPUT_VCF --hardy --keep ${LIST}/WP.list --out WP
awk 'NR>1' WP.hwe | awk '{split($3,el,"/"); if ((el[1]+el[2]+el[3])>0) print (el[2]/(el[1]+el[2]+el[3]))}' > WP.het.txt
vcftools --gzvcf $INPUT_VCF --hardy --keep ${LIST}/FU.list --out FU
awk 'NR>1' FU.hwe | awk '{split($3,el,"/"); if ((el[1]+el[2]+el[3])>0) print (el[2]/(el[1]+el[2]+el[3]))}' > FU.het.txt

# SFS
vcftools --gzvcf $INPUT_VCF --counts --keep ${LIST}/SP.list --out SP
bash ind.SFS.sh SP.frq.count > SP.SFS.txt
vcftools --gzvcf $INPUT_VCF --counts --keep ${LIST}/WP.list --out WP
bash ind.SFS.sh WP.frq.count > WP.SFS.txt
vcftools --gzvcf $INPUT_VCF --counts --keep ${LIST}/FU.list --out FU
bash ind.SFS.sh FU.frq.count > FU.SFS.txt

# TajimaD
vcftools --gzvcf $INPUT_VCF --TajimaD 10000 --keep ${LIST}/SP.list --out SP
awk 'NR>1' SP.Tajima.D | awk '{print $4}' > SP.TajimaD.txt
vcftools --gzvcf $INPUT_VCF --TajimaD 10000 --keep ${LIST}/WP.list --out WP
awk 'NR>1' WP.Tajima.D | awk '{print $4}' > WP.TajimaD.txt
vcftools --gzvcf $INPUT_VCF --TajimaD 10000 --keep ${LIST}/FU.list --out FU
awk 'NR>1' FU.Tajima.D | awk '{print $4}' > FU.TajimaD.txt
