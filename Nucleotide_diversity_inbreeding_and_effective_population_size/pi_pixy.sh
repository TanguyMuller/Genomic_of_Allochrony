#!/bin/bash
# Script to calcul Pi for each chromosome

# Extraire le nom du chromosome
chr=$(basename "$i" | cut -d'.' -f2)

# Créer un dossier par chromosome
mkdir -p /home/mullerta/work/SP_WP/Analyses/Nucleotide_diversity/"$chr"

# Indexer le fichier VCF
gvcf="Portugal.vcf.gz"
tabix -p vcf "${i}"

# Exécuter la commande pixy
pixy --stats pi \
   --vcf "${i}" \
   --populations /home/mullerta/work/SP_WP/Analyses/Nucleotide_diversity/population_list/population.portugal.list \
   --window_size 100000 \
   --n_cores 6 \
   --output_folder /home/mullerta/work/SP_WP/Analyses/Nucleotide_diversity/"$chr"

done
