#!/bin/bash 

# Script for Admixture analysis

vcf=$1
base_name=$(basename "${vcf}" .vcf.gz)
nbr_grp=$2
chr=$3

cd $chr

admixture --cv $base_name.bed $nbr_grp > log${nbr_grp}.out

### Command to run 7_Admixture.sh
#for i in {2..10}; do echo "sbatch 7_Admixture.sh ind.portugal.freebayes_filt_maf005_ldprune.vcf.gz $i Auto" ; done | bash
#for i in {2..10}; do echo "sbatch 7_Admixture.sh ind.portugal.freebayes_filt_maf005_ldprune.vcf.gz $i chrZ" ; done | bash
