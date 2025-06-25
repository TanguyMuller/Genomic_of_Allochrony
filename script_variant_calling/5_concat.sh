### Script for concat vcf file 

#! /bin/bash

# Note : Chunk were ordered by chromosome and position in all the lists. 
bcftools concat -f all.vcf.freebayes.list | bcftools view -S all.sample.ord - -Oz > all.portugal.freebayes.vcf.gz
bcftools concat -f pool.vcf.freebayes.list | bcftools view -S pool.sample.ord - -Oz > pool.portugal.freebayes.vcf.gz
bcftools concat -f ind.vcf.freebayes.list | bcftools view -S ind.sample.ord - -Oz > ind.portugal.freebayes.vcf.gz
bcftools concat -f all.vcf.freebayes.chrZ.list | bcftools view -S all.sample.ord - -Oz > all.portugal.freebayes.chrZ.vcf.gz
bcftools concat -f pool.vcf.freebayes.chrZ.list | bcftools view -S pool.sample.ord - -Oz > pool.portugal.freebayes.chrZ.vcf.gz
bcftools concat -f ind.vcf.freebayes.chrZ.list | bcftools view -S ind.sample.ord - -Oz > ind.portugal.freebayes.chrZ.vcf.gz
