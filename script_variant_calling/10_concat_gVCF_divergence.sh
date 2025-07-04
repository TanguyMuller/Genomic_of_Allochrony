#! /bin/bash
### Script to concat gvcf file divergence
# Note : gVCF files are ordered by chromosome in gvcf.divergence.portugal.list
bcftools concat -f gvcf.divergence.portugal.list | bcftools view -S sample.ord - -Oz > gVCF.divergence.vcf.gz
