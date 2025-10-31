### Script for create gVCF file of each chr 

#! /bin/bash

BAMLIST=$1
REGION=$2
PATH_TO_ASSEMBLY=/PATH/TO/ASSEMBLY

bcftools mpileup -Ou -C 50 --max-depth 500 -q 20 -Q 20 -a DP,AD -f $PATH_TO_ASSEMBLY -b $BAMLIST -r $REGION | bcftools call -m -Oz -o gVCF_Portugal/Portugal.${REGION}.gVCF.vcf.gz

### Command to run 7_gVCF_ind_portugal.sh
#awk '{print "bash 7_gVCF_ind_portugal.sh ind.bam.list ",$1}' chr.list | bash
