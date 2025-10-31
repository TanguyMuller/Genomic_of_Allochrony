### Script to concat gvcf file with the 2 most covered individuals of LSP, LWP, FU and outgroups

#! /bin/bash

# Note : gVCF files are ordered by chromosome in gvcf.divergence.portugal.list
bcftools concat -f gvcf.divergence.portugal.list | bcftools view -S sample.ord - -Oz > gVCF.divergence.vcf.gz
