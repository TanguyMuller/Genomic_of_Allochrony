#!/bin/bash

### Calcul of Pi using pixy
# Extraire le nom du chromosome
chr=$(basename "$i" | cut -d'.' -f2)
mkdir -p /home/mullerta/work/SP_WP/Analyses/Nucleotide_diversity/"$chr"
# Indexer le fichier VCF
tabix -p vcf "${i}"

# Exécuter la commande pixy
pixy --stats pi \
   --vcf "${i}" \
   --populations /home/mullerta/work/SP_WP/Analyses/Nucleotide_diversity/population_list/population.portugal.list \
   --window_size 100000 \
   --n_cores 6 \
   --output_folder /home/mullerta/work/SP_WP/Analyses/Nucleotide_diversity/"$chr"

done
