#!/bin/bash
# Annotate variants using SnpEff and index the resulting VCF

java -Xmx4g -jar snpEff/snpEff.jar PPM demo.portugal.vcf.gz | bgzip -c > load.ann.vcf.gz
tabix -p vcf load.ann.vcf.gz

java -Xmx4g -jar snpEff/snpEff.jar PPM Z.demo.portugal.vcf.gz | bgzip -c > Z_load.ann.vcf.gz
tabix -p vcf Z_load.ann.vcf.gz
