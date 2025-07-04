#! /bin/bash
### Script to concat gvcf file with portuguese individuals 
# Note : gVCF files are ordered by chromosome in gvcf.ind.portugal.list
bcftools concat -f gvcf.ind.portugal.list | bcftools view -S sample.ord - -Oz > gVCF.portugal.vcf.gz
