#!/bin/bash
# Script to calculate Dxy (absolute genetic divergence) between populations using pixy
# Processes VCF file and generates pairwise Dxy values in sliding windows

# Path to compressed gVCF file
gVCF="gVCF.divergence.vcf.gz"

# Index the gVCF file (required for pixy)
tabix -p vcf ${gVCF}

# Run pixy to calculate Dxy statistics
pixy --stats dxy \
	 --vcf ${gVCF} \
     --populations population.portugal.list \
     --window_size 100000 \
     --n_cores 6 \
     --output_folder calcul_dxy_by_pop

# Extract pairwise Dxy values for each population combination
# Population codes
values=("SP" "WP" "FU" "ES" "RI" "OK" "TU" "BO" "PI" "WI")

# Generate all pairwise combinations
for ((i = 1; i <= ${#values[@]}; i++)); do
  current_i=${values[i - 1]}
  
  for ((j = i + 1; j <= ${#values[@]}; j++)); do
    current_j=${values[j - 1]}
    
    # Create output file for each population pair
    output_file="calcul_dxy_by_pop/${current_i}_${current_j}_dxy.txt"
    
    # Add header to output file
    head -n 1 pixy_dxy.txt > "$output_file"
    
    # Extract Dxy values for current population pair
    awk -v pattern1="$current_i" -v pattern2="$current_j" '$1 ~ pattern1 && $2 == pattern2 {print $0}' calcul_dxy_by_pop/pixy_dxy.txt >> "$output_file"
  done
done
