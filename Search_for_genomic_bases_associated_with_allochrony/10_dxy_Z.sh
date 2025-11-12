#!/bin/bash
# Script for sliding window of Dxy on the Z chromosome

gVCF="Z.divergence.gVCF.vcf.gz"

# Indexer le fichier VCF
tabix -f -p vcf ${gVCF}

# Exécuter la commande pixy
pixy --stats dxy \
     --vcf ${gVCF} \
     --populations populations_divergence.txt \
     --window_size 2000 \
     --n_cores 6 \
     --output_folder Z/calcul_dyx_by_pop

values=("SP" "WP" "FU" "F1" "LSP")

for ((i = 1; i <= ${#values[@]}; i++)); do
  current_i=${values[i - 1]}
  
  for ((j = i + 1; j <= ${#values[@]}; j++)); do
    current_j=${values[j - 1]}
    output_file="calcul_dxy_by_pop_selection/${current_i}_${current_j}_dxy.txt"
    head -n 1 calcul_dxy_by_pop_selection/pixy_dxy.txt > "$output_file"
    awk -v pattern1="$current_i" -v pattern2="$current_j" '$1 ~ pattern1 && $2 == pattern2 {print $0}' calcul_dxy_by_pop_selection/pixy_dxy.txt >> "$output_file"
  done
done
