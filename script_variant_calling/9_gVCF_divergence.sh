#! /bin/bash

BAMLIST=$1
REGION=$2
PATH_TO_ASSEMBLY=/PATH/TO/ASSEMBLY

bcftools mpileup -Ou -C 50 --max-depth 500 -q 20 -Q 20 -a DP,AD -f $PATH_TO_ASSEMBLY -b $BAMLIST -r $REGION | bcftools call -m -Oz -o gVCF_Divergence/Divergence.${REGION}.gVCF.vcf.gz

### Command to run 9_gVCF_divergence.sh
#awk '{print "bash 9_gVCF_divergence.sh divergence.bam.list ",$1}' chr.list | bash
