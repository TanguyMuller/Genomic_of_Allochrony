### Script for concat vcf file 

#! /bin/bash

# Note : Chunk were ordered by chromosome and position in vcf.freebayes.list
bcftools concat -f vcf.freebayes.list | bcftools view -S sample.ord - -Oz > all.portugal.freebayes.vcf.gz
