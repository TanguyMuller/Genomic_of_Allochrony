#!/bin/bash
# Script to filter SNPs based on Minor Allele Frequency (MAF) threshold of 5%
# Filters out sites with MAF < 5% in all populations (keeps site if at least one population has MAF >= 5%)

VCF=$1
chr=$2
base_name=$(basename "$VCF" .vcf.gz)
path_list="/path/to/list"

# Initialize statistics file
echo "MAF 5% FILTERING STATISTICS FOR ${VCF}" >> filter.MAF005.snp.${chr}.txt
echo "----------------------------------------" >> filter.MAF005.snp.${chr}.txt
echo "Date: $(date)" >> filter.MAF005.snp.${chr}.txt
echo "" >> filter.MAF005.snp.${chr}.txt
echo "Number of SNPs in initial file: $(bcftools view -H $VCF | wc -l)" >> filter.MAF005.snp.${chr}.txt
echo "" >> filter.MAF005.snp.${chr}.txt
echo "Filtering sites with MAF<5% in all populations (if at least one population or pool has MAF>5%, keep the site)..." >> filter.MAF005.snp.${chr}.txt

# Define group sizes based on input file type
if [[ $base_name == "all.portugal.freebayes_filt" ]]; then
    group_size="25 18 10 4 4 4 4 4 9 8 1 1 1 1 1 1 1 1 1 1"  
elif [[ $base_name == "ind.portugal.freebayes_filt" ]]; then
    group_size="25 18 10 4 4 4 4 4 9 8"
else
    group_size="1 1 1 1 1 1 1 1 1 1"
fi

# Filter SNPs based on MAF threshold
zcat ${VCF} | awk -v group_size="$group_size" '
BEGIN { OFS = "\t" }
{
    if ($1 ~ /^#/) {
        # Skip header lines
        next
    }
    
    # Define group boundaries based on the number of individuals and pools
    split(group_size, sizes, " ")
    group_start = 10
    at_least_one_above_threshold = 0  # Flag to check if any ratio is >= 0.05
    at_least_one_below_threshold = 0  # Flag to check if any ratio is < 0.95
    
    for (g = 1; g <= length(sizes); g++) {
        group_end = group_start + sizes[g] - 1
        
        sum_dp = 0
        sum_ad_second = 0
        
        for (i = group_start; i <= group_end; i++) {
            split($i, fields, ":")  # Split fields by ":"
            dp = fields[2]           # DP is the second field
            ad = fields[3]           # AD is the third field
            split(ad, ad_values, ",")  # Split AD into individual allele depths
            sum_dp += dp
            sum_ad_second += ad_values[2]  # Sum the second value of AD for the group
        }
        
        # Calculate the AD/DP ratio for this group
        if (sum_dp > 0) {
            ratio = sum_ad_second / sum_dp
        } else {
            ratio = 0  # Handle cases where DP is 0, treating it as a ratio of 0
        }
        
        # Store the ratio for this group
        ratios[g] = ratio
        
        # Check if the ratio is greater than or equal to 0.05
        if (ratio >= 0.05) {
            at_least_one_above_threshold = 1
        }
        
        # Check if the ratio is less than 0.95
        if (ratio <= 0.95) {
            at_least_one_below_threshold = 1
        }
        
        group_start = group_end + 1  # Update the start for the next group
    }
    
    # Print only the chr:pos if both conditions are met
    if (at_least_one_above_threshold && at_least_one_below_threshold) {
        print $1 "\t" $2
    }
}
' > ${path_list}/${base_name}.positions.maf005.list

# Create filtered dataset
bcftools view -T ${path_list}/${base_name}.positions.maf005.list $VCF -Oz > ${base_name}_maf005.vcf.gz

echo "Number of SNPs after MAF<5% filtering: $(bcftools view -H ${base_name}_maf005.vcf.gz | wc -l)" >> filter.MAF005.snp.${chr}.txt
echo "" >> filter.MAF005.snp.${chr}.txt
echo "End filtering"

### Commands to run 1_filtering_maf005.sh
#bash 1_filtering_maf005.sh all.portugal.freebayes_filt.vcf.gz Auto
#bash 1_filtering_maf005.sh pool.portugal.freebayes_filt.vcf.gz Auto
#bash 1_filtering_maf005.sh ind.portugal.freebayes_filt.vcf.gz Auto
#bash 1_filtering_maf005.sh all.portugal.freebayes_filt.chrZ.vcf.gz chrZ
#bash 1_filtering_maf005.sh all.portugal.freebayes_filt.chrZ.vcf.gz chrZ
#bash 1_filtering_maf005.sh all.portugal.freebayes_filt.chrZ.vcf.gz chrZ
