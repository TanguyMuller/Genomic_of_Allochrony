#!/bin/bash
# Script for run FROH analyses for ROH>100kb and ROH>500Kb 

for size in "100" "500";do
	for pop in "LSP" "LWP" "FU" "VI" "VA" "CA" "TA" "GR";do
		mkdir -p ${size}kb
		VCF_INPUT="ind_masked_HWE.vcf.gz"
		OUTPUT="${size}kb/${pop}.ROH.${size}kb"
		
		./plink --allow-extra-chr \
				--chr-set 49 \
				--double-id \
				--homozyg \
				--homozyg-density 20 \
				--homozyg-gap 100 \
				--homozyg-kb $size \
				--homozyg-snp 50 \
				--homozyg-window-het 2 \
				--homozyg-window-missing 5 \
				--homozyg-window-snp 50 \
				--out $OUTPUT \
				--vcf $VCF_INPUT
	done;
done
