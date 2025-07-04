### Script to concat vcf file 

#! /bin/bash

# Note : gVCF files are ordered by chromosome in gvcf.ind.portugal.list
bcftools concat -f gvcf.ind.portugal.list | bcftools view -S sample.ord - -Oz > gVCF.portugal.vcf.gz
