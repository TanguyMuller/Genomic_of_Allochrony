### Script for concat vcf file 

#! /bin/bash

# Note : Chunk were ordered by chromosome and position in vcf.freebayes.list
bcftools concat -f all.vcf.freebayes.list | bcftools view -S all.sample.ord - -Oz > all.portugal.freebayes.vcf.gz
bcftools concat -f pool.vcf.freebayes.list | bcftools view -S pool.sample.ord - -Oz > pool.portugal.freebayes.vcf.gz
bcftools concat -f ind.vcf.freebayes.list | bcftools view -S ind.sample.ord - -Oz > ind.portugal.freebayes.vcf.gz
