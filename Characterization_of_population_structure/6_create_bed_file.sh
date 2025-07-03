#!/bin/bash 
# Script to generate BED files for Admixture analysis

vcf=$1
base_name=$(basename "${vcf}" .vcf.gz)
chr=$2

mkdir -p $chr
cd $chr

# Generate the input file in plink format
plink --vcf $vcf --make-bed --out $base_name --allow-extra-chr --chr-set 49 --double-id --keep-fam ind.portugal.list

# Change chr in 0 cause plink accept only human chr
awk '{$1="0";print $0}' $base_name.bim > $base_name.bim.tmp
mv $base_name.bim.tmp $base_name.bim

### Command
#sbatch 6_create_bed_file.sh ind.portugal.freebayes_filt_maf005_ldprune.vcf.gz Auto
#sbatch 6_create_bed_file.sh ind.portugal.freebayes.chrZ_filt_maf005_ldprune.vcf.gz chrZ 
