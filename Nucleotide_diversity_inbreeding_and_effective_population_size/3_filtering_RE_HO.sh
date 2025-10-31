#!/bin/bash 
# Script to filter SNPs showing high heterozygosity (Ho ≥ 70%)
# across all PPM individuals, or located within repeated regions
# identified in the reference genome.
#
# Command
# bash 3_filtering_RE_HO.sh ind.vcf.gz Auto

VCF_IND=$1
chr=$2
base_name=$(basename "$VCF_IND" .vcf.gz)
path_list="/your/way/to/list/"

# Initialiser le fichier de statistiques
echo "FILTERING STATISTICS: RE (repeats) and HO (heterozygosity) for ${VCF_IND}" > filter.RE.HO.snp.${chr}.txt
echo "----------------------------------------" >> filter.RE.HO.snp.${chr}.txt
echo "Date: $(date)" >> filter.RE.HO.snp.${chr}.txt
echo "" >> filter.RE.HO.snp.${chr}.txt
echo "Starting RE and HO filtering process..." >> filter.RE.HO.snp.${chr}.txt
echo "------------------------------------------------------------" >> filter.RE.HO.snp.${chr}.txt
echo "" >> filter.RE.HO.snp.${chr}.txt

echo "Calculating observed heterozygosity (Ho) per site..." >> filter.RE.HO.snp.${chr}.txt
echo "" >> filter.RE.HO.snp.${chr}.txt
vcftools --gzvcf $VCF_IND --keep ${path_list}/ind_list.list --hardy --out ${base_name}

# Remove sites located within repeated regions
echo "Filtering out SNPs located in transposable elements or other repeated regions..." >> filter.RE.HO.snp.${chr}.txt
bcftools view -T ^${path_list}/positions.RE.list ${base_name}.vcf.gz -Oz > ${base_name}_masked.vcf.gz
echo "Number of SNPs after RE filtering: $(gunzip -c ${base_name}_masked.vcf.gz | grep -c ^chr)" >> filter.RE.HO.snp.${chr}.txt
echo "" >> filter.RE.HO.snp.${chr}.txt

# Remove sites with >70% heterozygous individuals
echo "Filtering SNPs with more than 70% heterozygous individuals..." >> filter.RE.HO.snp.${chr}.txt
awk 'NR>1' ${base_name}.hwe | awk '{split($3,el,"/"); if ((el[1]+el[2]+el[3])>0 && (el[2]/(el[1]+el[2]+el[3]))>0.7) print $1"\t"$2}' > ${path_list}/${base_name}.positions.HWE.list
bcftools view -T ^${path_list}/${base_name}.positions.HWE.list ${base_name}_masked.vcf.gz -Oz > ${base_name}_masked_HWE.vcf.gz
bcftools index ${base_name}_masked_HWE.vcf.gz
echo "Number of SNPs after Ho>70% filtering: $(gunzip -c ${base_name}_masked_HWE.vcf.gz | grep -c ^chr)" >> filter.RE.HO.snp.${chr}.txt
rm ${base_name}.hwe ${base_name}_masked.vcf.gz

# Summary of filtering
echo "" >> filter.RE.HO.snp.${chr}.txt
echo "Final filtered file created: ${base_name}_masked_HWE.vcf.gz" >> filter.RE.HO.snp.${chr}.txt
echo "" >> filter.RE.HO.snp.${chr}.txt
echo "SUMMARY" >> filter.RE.HO.snp.${chr}.txt
echo "------" >> filter.RE.HO.snp.${chr}.txt
echo "1) Initial dataset (no filters) = ${base_name}.vcf.gz ($(gunzip -c ${base_name}.vcf.gz | grep -c ^chr) SNPs)" >> filter.RE.HO.snp.${chr}.txt
echo "2) Filtered dataset (RE + Ho): = ${base_name}_masked_HWE.vcf.gz ($(gunzip -c ${base_name}_masked_HWE.vcf.gz | grep -c ^chr) SNPs)" >> filter.RE.HO.snp.${chr}.txt

echo ""Filtering completed. See statistics in filter.RE.HO.snp.${chr}.txt"
