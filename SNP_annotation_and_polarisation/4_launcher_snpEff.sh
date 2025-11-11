#!/bin/bash
# Annotate variants using SnpEff and index the resulting VCF

java -Xmx4g -jar snpEff/snpEff.jar PPM demo.portugal.vcf.gz | bgzip -c > load.ann.vcf.gz
tabix -p vcf load.ann.vcf.gz
