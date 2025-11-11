#!/bin/bash
# Extract the fields for next load analyses

# Parameters
VCF=$1
OUT=$2

# Run
#freebayes:
java -jar snpEff/SnpSift.jar extractFields -e "NA" vcf/${VCF} "CHROM" "POS" "REF" "ALT" "GEN[*].GT" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].HGVS_P" | sed 's/,/\t/g' > vcf/${OUT}

#Command
#sbatch 1_create_table_GT.sh load.ann.vcf.gz load.ann.GT.txt
#sbatch 1_create_table_GT.sh Z.load.ann.vcf.gz Z.load.ann.GT.txt
