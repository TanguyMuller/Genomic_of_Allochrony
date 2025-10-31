Nouvelle conversation

j'ai fais library(leastcostpath)
library(terra)
library(sf)
library(rnaturalearth)
library(elevatr)
library(dplyr)

pops <- data.frame(
POPULATIONS = c("ANG","AVI","BAR","BAY","BEA","BIL","BOU","BRE","CER",
"CES","CHA","CTB","CTH","ERQ","ETA","FIR","FRO","LAR","LLE",
"LPR","LYO","MAI","MTP","ORL","PAM","PAR","PIE","PLO","ROD",
"ROL","TAR","THO","TOU","TRE","TRO","VES","VOU","VTB","VTH"),
LATITUDE = c(47.324,44.155,41.573,43.411,47.027,43.243,46.916,48.236,42.449,
44.749,48.897,37.085,37.087,48.638,48.440,43.944,42.502,45.130,41.663,
46.254,45.654,48.451,43.631,47.790,42.511,44.610,45.150,47.568,44.480,
43.945,40.890,46.350,43.762,45.518,48.297,45.618,46.600,44.255,44.355),
LONGITUDE = c(-0.609,5.146,2.389,-1.351,4.807,-2.718,2.300,-4.482,3.142,
-0.788,4.088,-3.470,-3.451,-2.446,2.286,10.912,2.099,1.257,0.639,
-1.516,4.786,1.445,3.839,1.724,-1.674,10.072,7.053,-3.134,2.464,
8.129,0.679,6.356,1.360,11.167,4.074,11.177,0.170,5.196,5.210)
)
#
#
#

# Initialiser la carte
plot(cost_surface, main = "Least Cost Paths depuis CTB")

# Initialiser un data.frame pour stocker les résultats
results <- data.frame(Population = character(), Distance_km = numeric(), stringsAsFactors = FALSE)

# Boucle sur chaque population (sauf CTB)
for (i in 1:nrow(pops)) {
if (pops$POPULATIONS[i] != "CTB") {
# Coordonnées de CTB (origine)
coords_ctb <- data.frame(
ID = "CTB",
LAT = pops$LATITUDE[pops$POPULATIONS == "CTB"],
LON = pops$LONGITUDE[pops$POPULATIONS == "CTB"]
)

# Coordonnées de la population cible
coords_target <- data.frame(
ID = pops$POPULATIONS[i],
LAT = pops$LATITUDE[i],
LON = pops$LONGITUDE[i]
)

# Combine les deux points
points_df <- rbind(coords_ctb, coords_target)

# Crée les points SF
points_sf <- st_as_sf(points_df, coords = c("LON", "LAT"), crs = 4326)
points_proj <- st_transform(points_sf, 3035)

# Calculer le LCP
cs <- create_cs(cost_surface)
lcps <- create_lcp(x = cs, origin = points_proj[1, ], destination = points_proj[2, ], cost_distance = TRUE)

# Calculer la distance
dist_km <- as.numeric(st_length(lcps)) / 1000

# Stocker le résultat
results <- rbind(results, data.frame(Population = pops$POPULATIONS[i], Distance_km = round(dist_km, 2)))

# Visualisation
plot(st_geometry(lcps), add = TRUE, col = "red", lwd = 1)
plot(st_geometry(points_proj), add = TRUE, col = "magenta", pch = 19, cex=0.2)
#points_sf <- st_as_sf(coords, coords=c("LON","LAT"), crs=4326)
}
}

# Afficher le tableau des distances
print(results)

# Exporter le tableau
write.csv(results, "distances_lcp.csv", row.names = FALSE)
et je voudrais ecrire un materiel et méthodes en anglais de ce que j'ai fait en assez concis

Here's a concise Materials and Methods section in English for what you've done:
Materials and Methods

Least Cost Path Analysis

We conducted a least cost path (LCP) analysis to estimate travel distances between the CTB population (Cantabria, Spain; 37.085°N, 3.470°W) and 38 other European populations.

Digital elevation data were obtained using the elevatr package in R (v. X.X). A cost surface was generated from the elevation raster, where terrain difficulty influences travel cost. All LCP calculations were performed using the leastcostpath package (v. X.X) in R (R Core Team, 2024).

Population coordinates were transformed from WGS84 (EPSG:4326) to ETRS89/LAEA Europe projection (EPSG:3035) using the sf package. For each population pair, we computed the least cost path from CTB to the target population using the create_cs() and create_lcp() functions with cost-distance weighting enabled.

Path distances were calculated as the geodesic length of the resulting polyline geometries and converted to kilometers. All spatial analyses were conducted in R version X.X.X (R Core Team, 2024) using packages terra, sf, dplyr, and rnaturalearth for data manipulation and visualization.

Notes:

    Replace "v. X.X" with actual package versions
    Add more detail about your cost surface parameterization if you used specific cost functions (slope, Tobler's hiking function, etc.)
    Consider adding a sentence about why you chose CTB as the origin point if relevant to your study

ah mais beaucoup plus concis t'es pas obligé de mettre toute les étapes c'est du global qu'il faut faire

Voici une version beaucoup plus concise :
Materials and Methods

Least Cost Path Analysis

Least cost path (LCP) distances were calculated between the CTB population (Cantabria, Spain) and 38 other European populations using digital elevation data. Cost surfaces were generated from elevation rasters, and optimal paths were computed using the leastcostpath package in R. All coordinates were projected to ETRS89/LAEA Europe (EPSG:3035), and distances were measured in kilometers along the resulting paths.

Ou encore plus court si besoin :

Least Cost Path Analysis

Least cost path distances between CTB and 38 European populations were calculated using digital elevation data and the leastcostpath R package. Paths were computed in ETRS89/LAEA Europe projection (EPSG:3035) and distances measured in kilometers.

We conducted a least cost path (LCP) analysis to estimate geographical distances between the CTB population and 38 other European populations.
Digital elevation data were obtained using the \textit{elevatr} package in R (v. X.X). A cost surface was generated from the elevation raster, where terrain difficulty influences travel cost. All LCP calculations were performed using the \textit{leastcostpath} package (v. X.X) in R. For each population pair, we computed the least cost path from CTB to the target population using the \textit{create\_cs} and \textit{create\_lcp} functions with cost-distance weighting enabled. Path distances were calculated as the geodesic length of the resulting polyline geometries and converted to kilometers.  ça s'est parfait juste il ya trop de phrase tu peux merger certaines surtout au debut

Voici une version plus fluide avec les phrases du début fusionnées :

We conducted a least cost path (LCP) analysis to estimate geographical distances between the CTB population and 38 other European populations using digital elevation data obtained with the \textit{elevatr} package in R (v. X.X). A cost surface was generated from the elevation raster, where terrain difficulty influences travel cost, and all LCP calculations were performed using the \textit{leastcostpath} package (v. X.X) in R. For each population pair, we computed the least cost path from CTB to the target population using the \textit{create_cs} and \textit{create_lcp} functions with cost-distance weighting enabled. Path distances were calculated as the geodesic length of the resulting polyline geometries and converted to kilometers.

ok autre chose j'ai un code qui ressemble à ça module purge
module load bioinfo/mosdepth/0.3.6
mosdepth -t 12 -nx $BAMdir/stats/${rg} $BAMdir/${rg}.bam
rm $BAMdir/stats/${rg}.mosdepth.global.dist.txt
 et je veux enfait lancer mosdepth sur des bam dans un fichier txt (une liste des bam) et je veux donner en resultats des noms qui sont dans un autre fichier

Voici une solution avec une boucle qui lit les deux fichiers en parallèle :
bash

module purge
module load bioinfo/mosdepth/0.3.6

# Lire les fichiers ligne par ligne
paste bam_list.txt names_list.txt | while read bam_file output_name; do
    mosdepth -t 12 -nx $BAMdir/stats/${output_name} ${bam_file}
    rm $BAMdir/stats/${output_name}.mosdepth.global.dist.txt
done

Ou si tu préfères avec un index :
bash

module purge
module load bioinfo/mosdepth/0.3.6

# Lire les fichiers dans des arrays
mapfile -t bam_files < bam_list.txt
mapfile -t output_names < names_list.txt

# Boucle sur les indices
for i in "${!bam_files[@]}"; do
    mosdepth -t 12 -nx $BAMdir/stats/${output_names[$i]} ${bam_files[$i]}
    rm $BAMdir/stats/${output_names[$i]}.mosdepth.global.dist.txt
done

Notes :

    Les deux fichiers doivent avoir le même nombre de lignes
    La première solution avec paste est plus robuste pour les noms avec espaces
    Remplace bam_list.txt et names_list.txt par tes fichiers

sur plusieurs coeurs si c'est possible car en boucle ça va prendre du temps

Voici une solution avec GNU Parallel pour paralléliser sur plusieurs cœurs :
bash

module purge
module load bioinfo/mosdepth/0.3.6

# Avec GNU Parallel (recommandé)
parallel -j 4 --link 'mosdepth -t 12 -nx $BAMdir/stats/{2} {1} && rm $BAMdir/stats/{2}.mosdepth.global.dist.txt' :::: bam_list.txt names_list.txt

Ou si GNU Parallel n'est pas disponible, avec xargs :
bash

module purge
module load bioinfo/mosdepth/0.3.6

paste bam_list.txt names_list.txt | xargs -P 4 -n 2 bash -c 'mosdepth -t 12 -nx $BAMdir/stats/$1 $0 && rm $BAMdir/stats/$1.mosdepth.global.dist.txt'

Ou dans un script SLURM avec job array :
bash

#!/bin/bash
#SBATCH --array=1-N  # N = nombre de BAMs
#SBATCH --cpus-per-task=12

module purge
module load bioinfo/mosdepth/0.3.6

bam_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bam_list.txt)
output_name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" names_list.txt)

mosdepth -t 12 -nx $BAMdir/stats/${output_name} ${bam_file}
rm $BAMdir/stats/${output_name}.mosdepth.global.dist.txt

Notes :

    -j 4 ou -P 4 : nombre de jobs en parallèle (ajuste selon tes ressources)
    Attention : chaque mosdepth utilise déjà 12 threads (-t 12), donc si tu lances 4 en parallèle = 48 threads totaux
    La solution SLURM job array est la meilleure pour un cluster

ok mtn je veux aller dans ce dossier tu prend tous les fichiers .mosdepth.summary.txt et en premier tu fais head -50 ensuite tu prends la 4 colonne et tu marques dans un fichier DP.txt en premier ce qui a avant .mosdepth.summary.txt et enusite tu met la moyenne des ligne 1 à 49 et dans la derniere colonne tu met la 50ème colonne

Voici un script bash pour faire ça :
bash

#!/bin/bash

# Créer/vider le fichier de sortie
> DP.txt

# Boucle sur tous les fichiers .mosdepth.summary.txt
for file in *.mosdepth.summary.txt; do
    # Extraire le nom (sans .mosdepth.summary.txt)
    name=$(basename "$file" .mosdepth.summary.txt)
    
    # Prendre les 50 premières lignes, extraire colonne 4
    col4=$(head -50 "$file" | awk '{print $4}')
    
    # Calculer la moyenne des lignes 1-49 et extraire la ligne 50
    mean=$(echo "$col4" | head -49 | awk '{sum+=$1} END {print sum/NR}')
    line50=$(echo "$col4" | sed -n '50p')
    
    # Écrire dans le fichier de sortie
    echo -e "${name}\t${mean}\t${line50}" >> DP.txt
done

echo "Fichier DP.txt créé avec succès!"

Ou version plus compacte :
bash

for file in *.mosdepth.summary.txt; do
    name=${file%.mosdepth.summary.txt}
    head -50 "$file" | awk -v name="$name" 'NR<=49 {sum+=$4; count++} NR==50 {print name"\t"sum/count"\t"$4}' >> DP.txt
done

Le fichier DP.txt aura 3 colonnes :

    Nom du sample
    Moyenne de la colonne 4 (lignes 1-49)
    Valeur de la colonne 4 (ligne 50)

juste pardon c'est pas head -50 mais les lignes 2 à 51

Voici la version corrigée pour les lignes 2 à 51 :
bash

#!/bin/bash

# Créer/vider le fichier de sortie
> DP.txt

# Boucle sur tous les fichiers .mosdepth.summary.txt
for file in *.mosdepth.summary.txt; do
    # Extraire le nom (sans .mosdepth.summary.txt)
    name=$(basename "$file" .mosdepth.summary.txt)
    
    # Prendre les lignes 2 à 51, extraire colonne 4
    col4=$(sed -n '2,51p' "$file" | awk '{print $4}')
    
    # Calculer la moyenne des lignes 2-50 et extraire la ligne 51
    mean=$(echo "$col4" | head -49 | awk '{sum+=$1} END {print sum/NR}')
    line51=$(echo "$col4" | tail -1)
    
    # Écrire dans le fichier de sortie
    echo -e "${name}\t${mean}\t${line51}" >> DP.txt
done

echo "Fichier DP.txt créé avec succès!"

Ou version plus compacte :
bash

for file in *.mosdepth.summary.txt; do
    name=${file%.mosdepth.summary.txt}
    sed -n '2,51p' "$file" | awk -v name="$name" 'NR<=49 {sum+=$4; count++} NR==50 {print name"\t"sum/count"\t"$4}' >> DP.txt
done

Le fichier DP.txt aura 3 colonnes :

    Nom du sample
    Moyenne de la colonne 4 (lignes 2-50)
    Valeur de la colonne 4 (ligne 51)

je veux faire un MSMC2 et je suis à l'étape ou je dois préparer mes vcf et mes .bed mais je me demande si je dois prendre des individus de mon VCF ou mon gVCF

Pour MSMC2, tu dois utiliser ton VCF (pas le gVCF), mais avec quelques nuances importantes :

Pourquoi le VCF et pas le gVCF ?

    Le gVCF contient tous les sites (variants + non-variants) ce qui est énorme et non nécessaire
    MSMC2 a besoin des sites variants pour l'analyse + un masque de couverture (fichiers .bed) pour savoir quelles régions sont fiables

Ce dont MSMC2 a besoin :

    VCF phasé avec les variants de tes individus
        Doit être phasé (avec beagle, shapeit, whatshapp, etc.)
        Format : un VCF par chromosome
    Fichiers .bed de masque
        Ces fichiers indiquent les régions appelables (bonne couverture, pas répétées, etc.)
        C'est ici que ton gVCF peut être utile : tu peux l'utiliser pour générer le masque en identifiant les régions avec bonne couverture et qualité

Workflow typique :
bash

# 1. Utilise ton VCF phasé pour les variants
# 2. Génère le masque .bed depuis ton BAM ou gVCF (régions avec DP suffisant, qualité OK)
# 3. Lance MSMC2 avec VCF + masque

Tu es à quelle étape exactement ? Tu as déjà phasé ton VCF ?

ok j'ai un vcf contenant plusieurs individu et plusieurs chr et qui est dejà filtrer sur la couverture avec chaque individu qui a une couverture d'au moins 5  mais du coup c'est des filtres sur les sites variants

Ah d'accord, donc tu as un VCF multi-échantillons déjà filtré. Voici ce que tu dois faire pour MSMC2 :

Le problème actuel :

    Ton VCF est filtré sur les sites variants (DP ≥ 5)
    Mais MSMC2 a aussi besoin de savoir où sont les régions appelables (même sans variants)
    C'est pour ça qu'il faut les fichiers .bed de masque

Ce qu'il te faut faire :

    Séparer ton VCF par chromosome (si pas déjà fait)

bash

for chr in {1..22}; do
    bcftools view -r chr${chr} input.vcf.gz -O z -o chr${chr}.vcf.gz
done

    Phaser ton VCF (obligatoire pour MSMC2)

bash

# Avec beagle par exemple
java -jar beagle.jar gt=chr${chr}.vcf.gz out=chr${chr}.phased

    Générer les masques .bed depuis tes BAMs

bash

# Pour chaque individu, créer le masque des régions avec DP ≥ 5
# Tu peux utiliser mosdepth ou bcftools pour ça
samtools depth -a sample.bam | awk '$3 >= 5 {print $1"\t"$2-1"\t"$2}' | bedtools merge > sample.mask.bed

    Lancer MSMC2 avec VCF phasé + masques

Tu as déjà tes BAMs disponibles pour générer les masques ?

il faut pas aussi par individus ?? et du coup c'est pas grave si il y a des sites monorphes quand tu filtres sur 1 individus ??

Oui exactement, tu dois séparer par individu !

Pour MSMC2 il te faut :

    VCF par individu ET par chromosome

bash

for ind in sample1 sample2 sample3; do
    for chr in {1..22}; do
        bcftools view -s ${ind} -r chr${chr} input.vcf.gz -O z -o ${ind}.chr${chr}.vcf.gz
    done
done

    Sites monomorphes : c'est NORMAL et OK !
        Quand tu extrais un seul individu, beaucoup de sites deviennent monomorphes (0/0 ou 1/1)
        MSMC2 a besoin de ces sites monomorphes pour calculer correctement les temps de coalescence
        Ne les filtre surtout pas !
        Le masque .bed sert justement à dire "dans ces régions, l'absence de variant est informative"

Workflow complet :
bash

# Par individu
for ind in sample1 sample2; do
    # 1. Extraire l'individu (garde tous les sites, même monomorphes)
    bcftools view -s ${ind} -r chr1 multisample.vcf.gz -O z -o ${ind}.chr1.vcf.gz
    
    # 2. Phaser (peut skipper si déjà phasé)
    beagle gt=${ind}.chr1.vcf.gz out=${ind}.chr1.phased
    
    # 3. Masque depuis BAM
    samtools depth -a ${ind}.bam -r chr1 | awk '$3 >= 5' | bedtools merge > ${ind}.chr1.mask.bed
done

Tu veux faire du MSMC2 en within-population (un individu) ou cross-population (paires) ?

il fait la différence entre les sites 0/0 et les sites qui sont callable mais pas dans le VCF mais qui doit interpreter 0/0 du coup ??

Excellente question ! C'est exactement le rôle du masque .bed :

MSMC2 distingue 3 types de sites :

    Sites dans le VCF avec 0/0 → Région callable, génotype homozygote référence
    Sites dans le VCF avec 0/1 ou 1/1 → Région callable avec variant
    Sites PAS dans le VCF MAIS dans le masque .bed → Région callable, interprété comme 0/0
    Sites ni dans VCF ni dans masque → Région NON callable (ignorée)

Le masque .bed dit à MSMC2 :

    "Dans ces régions, j'ai une bonne couverture"
    "Si un site n'est pas dans le VCF mais est dans le masque → c'est un 0/0"
    "Si un site n'est ni dans le VCF ni dans le masque → ignore-le (mauvaise qualité)"

Donc en pratique :

    Ton VCF peut ne contenir QUE les variants (pas les 0/0) → plus léger
    Le masque .bed comble les trous et dit "ici c'était callable même sans variant"
    MSMC2 reconstruit : variant dans VCF = hétéro/homo alt, pas de variant mais masque = homo ref

C'est pour ça qu'il est CRUCIAL d'avoir un bon masque .bed !

Ton VCF actuel contient les sites 0/0 ou seulement les variants ?

et pour les 1/1 car c'est un ratio de hetero / homo donc 0/0 ou 1/1 c'est pareil ??

Excellente remarque ! Non, ce n'est PAS pareil pour MSMC2 !

MSMC2 utilise la densité d'hétérozygotes le long du génome :

    0/1 (hétérozygote) → Indique une coalescence récente entre les deux haplotypes
    0/0 et 1/1 (homozygotes) → Indique que les deux haplotypes ont le même allèle (coalescence plus ancienne ou récente selon le contexte)

Mais attention :

Pour MSMC2, ce qui compte c'est le pattern des hétérozygotes sur les haplotypes phasés. Après phasage, ton VCF devient :

    0|0 ou 1|1 → Les deux chromosomes sont identiques à ce site
    0|1 ou 1|0 → Les deux chromosomes diffèrent

MSMC2 regarde la distribution spatiale des différences entre haplotypes, pas vraiment si c'est 0/0 vs 1/1.

En pratique :

    Un site 0/0 dans ton VCF individuel = les deux chromosomes ont l'allèle référence
    Un site 1/1 dans ton VCF individuel = les deux chromosomes ont l'allèle alternatif
    Pour MSMC2, les deux sont "concordants" (pas de différence entre haplotypes)
    Seuls les 0/1 (ou 0|1, 1|0 après phasage) sont informatifs pour les coalescences récentes

Donc oui, 0/0 et 1/1 sont traités similairement (homozygotes), mais il faut quand même les garder dans le VCF !

oui donc je peux filtrer les 0/0 et les 1/1 tant que je laisse les 0/1

NON, ne filtre pas les 1/1 !

Voici pourquoi :

Les 1/1 sont INFORMATIFS pour MSMC2 :

    Un site 1/1 signifie que ton individu a l'allèle alternatif sur les deux chromosomes
    Cela donne de l'information sur la fréquence allélique et la structure haplotypique locale
    MSMC2 utilise le pattern complet des génotypes le long du génome, pas juste les hétérozygotes isolés

Tu PEUX filtrer les 0/0 (sites monomorphes référence) parce que :

    Ils sont implicites : si un site est dans le masque .bed mais pas dans le VCF → c'est un 0/0
    Ça allège le VCF

Donc le bon filtrage :
bash

# Garder UNIQUEMENT 0/1 et 1/1 (enlever les 0/0)
bcftools view -i 'GT="0/1" || GT="1/1"' input.vcf.gz -O z -o filtered.vcf.gz

Ou plus simple, garde tous les sites variants (AC > 0) :
bash

bcftools view -e 'AC==0' input.vcf.gz -O z -o variants_only.vcf.gz

Résumé : Garde 0/1 ET 1/1, retire seulement les 0/0 (qui seront implicites via le masque)

ok donc j'ai 4 individus qui sont CTB_1 BRE_1 BEA_1 et FIR_1 et je voudrais faire ça juste sur le chrZ

Parfait ! Voici le workflow complet pour MSMC2 sur le chrZ avec tes 4 individus :
bash

#!/bin/bash

# Variables
VCF="input.vcf.gz"  # Ton VCF multi-échantillons
CHR="chrZ"
SAMPLES=("CTB_1" "BRE_1" "BEA_1" "FIR_1")

# 1. Extraire chaque individu pour chrZ (garder variants seulement)
for sample in "${SAMPLES[@]}"; do
    bcftools view -s ${sample} -r ${CHR} ${VCF} | \
    bcftools view -e 'GT="0/0"' -O z -o ${sample}.${CHR}.vcf.gz
    bcftools index ${sample}.${CHR}.vcf.gz
done

# 2. Phaser avec Beagle (si pas déjà phasé)
for sample in "${SAMPLES[@]}"; do
    java -jar beagle.jar gt=${sample}.${CHR}.vcf.gz out=${sample}.${CHR}.phased
done

# 3. Générer les masques .bed depuis les BAMs
for sample in "${SAMPLES[@]}"; do
    # Ajuste le path vers ton BAM
    samtools depth -a -r ${CHR} ${sample}.bam | \
    awk '$3 >= 5 {print $1"\t"$2-1"\t"$2}' | \
    bedtools merge > ${sample}.${CHR}.mask.bed
done

# 4. Préparer pour MSMC2 (génère l'input au format MSMC)
# Tu devras utiliser generate_multihetsep.py fourni avec MSMC2

Questions :

    Tu as les fichiers BAM pour générer les masques ?
    Ton VCF est déjà phasé ou pas ?
    Tu veux faire du MSMC2 individuel (chaque individu séparément) ou cross-population (paires d'individus) ?

#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH -p workq
#SBATCH --mem=24G
module purge
module load bioinfo/Bcftools/1.21
module load bioinfo/Beagle/5.4
module load bioinfo/samtools/1.21
# Variables
VCF="freebayes_filt_0pct_missing_dp_filtered_masked_HWE.vcf.gz"
CHR="chr1"
SAMPLES=("CTB_1" "BRE_1" "BEA_1" "FIR_1")
BAM=("/work/project/loadexp/full_loadexp_data/bam/ES_Cor_i12" "/work/project/loadexp/full_loadexp_data/bam/FR_Bre_i01" "/work/project/loadexp/full_loadexp_data/bam/FR_Bea_i01" "/work/project/loadexp/full_loadexp_data/bam/IT_Fir_i01")
# 1. Extraire chaque individu pour chrZ (garder variants seulement)
for sample in "${SAMPLES[@]}"; do
    bcftools view -s ${sample} -r ${CHR} ${VCF} | \
    bcftools view -e 'GT="0/0"' -O z -o ${sample}.${CHR}.vcf.gz
    bcftools index ${sample}.${CHR}.vcf.gz
done
# 2. Phaser avec Beagle (si pas déjà phasé)
for sample in "${SAMPLES[@]}"; do
    java -jar /usr/local/bioinfo/src/Beagle/Beagle_v5.4/beagle.22Jul22.46e.jar gt=${sample}.${CHR}.vcf.gz out=${sample}.${CHR}.phased
done
# 3. Générer les masques .bed depuis les BAMs
for sample in "${BAM[@]}"; do
    # Ajuste le path vers ton BAM
    samtools depth -a -r ${CHR} ${sample}.bam | \
    awk '$3 >= 5 {print $1"\t"$2-1"\t"$2}' | \
    bedtools merge > ${SAMPLE}.${CHR}.mask.bed
done
 ça c'est le code que je veux faire pardon

Voici ton code corrigé :
bash

#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH -p workq
#SBATCH --mem=24G
module purge
module load bioinfo/Bcftools/1.21
module load bioinfo/Beagle/5.4
module load bioinfo/samtools/1.21
module load bioinfo/BEDTools/2.31.1  # Ajouter bedtools

# Variables
VCF="freebayes_filt_0pct_missing_dp_filtered_masked_HWE.vcf.gz"
CHR="chr1"
SAMPLES=("CTB_1" "BRE_1" "BEA_1" "FIR_1")
BAMS=("/work/project/loadexp/full_loadexp_data/bam/ES_Cor_i12.bam" \
      "/work/project/loadexp/full_loadexp_data/bam/FR_Bre_i01.bam" \
      "/work/project/loadexp/full_loadexp_data/bam/FR_Bea_i01.bam" \
      "/work/project/loadexp/full_loadexp_data/bam/IT_Fir_i01.bam")

# 1. Extraire chaque individu pour chr1 (garder variants seulement)
for sample in "${SAMPLES[@]}"; do
    bcftools view -s ${sample} -r ${CHR} ${VCF} | \
    bcftools view -e 'GT="0/0"' -O z -o ${sample}.${CHR}.vcf.gz
    bcftools index ${sample}.${CHR}.vcf.gz
done

# 2. Phaser avec Beagle
for sample in "${SAMPLES[@]}"; do
    java -jar /usr/local/bioinfo/src/Beagle/Beagle_v5.4/beagle.22Jul22.46e.jar \
    gt=${sample}.${CHR}.vcf.gz out=${sample}.${CHR}.phased
done

# 3. Générer les masques .bed depuis les BAMs
for i in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$i]}"
    bam="${BAMS[$i]}"
    
    samtools depth -a -r ${CHR} ${bam} | \
    awk '$3 >= 5 {print $1"\t"$2-1"\t"$2}' | \
    bedtools merge > ${sample}.${CHR}.mask.bed
done

echo "Done!"

Corrections apportées :

    Ajout du module bedtools
    Extension .bam ajoutée aux paths BAM
    Correction de la boucle BAM : utilisation de "${!SAMPLES[@]}" pour indexer correctement
    Correction ${SAMPLE} → ${sample} (majuscule/minuscule)

ok niquel mais juste j'au une liste sous forme chr\tposition pour les régions répétés que je veux enlever de mon bed car nn callable je fais comment

Tu peux soustraire ces régions de ton masque .bed avec bedtools subtract :
bash

#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH -p workq
#SBATCH --mem=24G
module purge
module load bioinfo/Bcftools/1.21
module load bioinfo/Beagle/5.4
module load bioinfo/samtools/1.21
module load bioinfo/BEDTools/2.31.1

# Variables
VCF="freebayes_filt_0pct_missing_dp_filtered_masked_HWE.vcf.gz"
CHR="chr1"
REPEATS="repeats.txt"  # Ton fichier chr\tposition
SAMPLES=("CTB_1" "BRE_1" "BEA_1" "FIR_1")
BAMS=("/work/project/loadexp/full_loadexp_data/bam/ES_Cor_i12.bam" \
      "/work/project/loadexp/full_loadexp_data/bam/FR_Bre_i01.bam" \
      "/work/project/loadexp/full_loadexp_data/bam/FR_Bea_i01.bam" \
      "/work/project/loadexp/full_loadexp_data/bam/IT_Fir_i01.bam")

# Convertir repeats en format BED (chr start end)
awk '{print $1"\t"$2-1"\t"$2}' ${REPEATS} | sort -k1,1 -k2,2n > repeats.bed

# 1. Extraire chaque individu pour chr1 (garder variants seulement)
for sample in "${SAMPLES[@]}"; do
    bcftools view -s ${sample} -r ${CHR} ${VCF} | \
    bcftools view -e 'GT="0/0"' -O z -o ${sample}.${CHR}.vcf.gz
    bcftools index ${sample}.${CHR}.vcf.gz
done

# 2. Phaser avec Beagle
for sample in "${SAMPLES[@]}"; do
    java -jar /usr/local/bioinfo/src/Beagle/Beagle_v5.4/beagle.22Jul22.46e.jar \
    gt=${sample}.${CHR}.vcf.gz out=${sample}.${CHR}.phased
done

# 3. Générer les masques .bed depuis les BAMs et enlever les répétitions
for i in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$i]}"
    bam="${BAMS[$i]}"
    
    samtools depth -a -r ${CHR} ${bam} | \
    awk '$3 >= 5 {print $1"\t"$2-1"\t"$2}' | \
    bedtools merge | \
    bedtools subtract -a stdin -b repeats.bed > ${sample}.${CHR}.mask.bed
done

echo "Done!"

Ce qui a été ajouté :

    Conversion de ton fichier repeats en format BED
    bedtools subtract -a stdin -b repeats.bed pour enlever les régions répétées du masque

parfait dernier truc juste je veux aussi que le chrZ soit pas callable

Voici la version avec chrZ exclu :
bash

#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH -p workq
#SBATCH --mem=24G
module purge
module load bioinfo/Bcftools/1.21
module load bioinfo/Beagle/5.4
module load bioinfo/samtools/1.21
module load bioinfo/BEDTools/2.31.1

# Variables
VCF="freebayes_filt_0pct_missing_dp_filtered_masked_HWE.vcf.gz"
CHR="chr1"
REPEATS="repeats.txt"  # Ton fichier chr\tposition
SAMPLES=("CTB_1" "BRE_1" "BEA_1" "FIR_1")
BAMS=("/work/project/loadexp/full_loadexp_data/bam/ES_Cor_i12.bam" \
      "/work/project/loadexp/full_loadexp_data/bam/FR_Bre_i01.bam" \
      "/work/project/loadexp/full_loadexp_data/bam/FR_Bea_i01.bam" \
      "/work/project/loadexp/full_loadexp_data/bam/IT_Fir_i01.bam")

# Convertir repeats en format BED (chr start end)
awk '{print $1"\t"$2-1"\t"$2}' ${REPEATS} > repeats.bed

# Créer un fichier BED pour chrZ (tout le chromosome)
# Ajuste la taille selon ton génome de référence
echo -e "chrZ\t0\t999999999" > chrZ.bed

# Combiner repeats + chrZ
cat repeats.bed chrZ.bed | sort -k1,1 -k2,2n > regions_to_exclude.bed

# 1. Extraire chaque individu pour chr1 (garder variants seulement)
for sample in "${SAMPLES[@]}"; do
    bcftools view -s ${sample} -r ${CHR} ${VCF} | \
    bcftools view -e 'GT="0/0"' -O z -o ${sample}.${CHR}.vcf.gz
    bcftools index ${sample}.${CHR}.vcf.gz
done

# 2. Phaser avec Beagle
for sample in "${SAMPLES[@]}"; do
    java -jar /usr/local/bioinfo/src/Beagle/Beagle_v5.4/beagle.22Jul22.46e.jar \
    gt=${sample}.${CHR}.vcf.gz out=${sample}.${CHR}.phased
done

# 3. Générer les masques .bed depuis les BAMs et enlever les répétitions + chrZ
for i in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$i]}"
    bam="${BAMS[$i]}"
    
    samtools depth -a -r ${CHR} ${bam} | \
    awk '$3 >= 5 {print $1"\t"$2-1"\t"$2}' | \
    bedtools merge | \
    bedtools subtract -a stdin -b regions_to_exclude.bed > ${sample}.${CHR}.mask.bed
done

echo "Done!"

Ajouts :

    Création d'un fichier chrZ.bed avec tout le chrZ
    Fusion avec repeats.bed dans regions_to_exclude.bed
    Soustraction des deux types de régions en une seule commande

mais dans le repeats il y a deja des sites dans le chrZ ça va pas créer de probleme de doublon

Non, pas de problème ! bedtools subtract gère automatiquement les chevauchements et doublons. Même si chrZ apparaît deux fois (une fois en entier et des positions spécifiques dans repeats), il va correctement enlever toutes les régions.

Mais si tu veux être plus propre, tu peux faire comme ça :
bash

# Convertir repeats en format BED
awk '{print $1"\t"$2-1"\t"$2}' ${REPEATS} > repeats.bed

# Ajouter chrZ entier seulement s'il n'y a pas déjà de lignes chrZ dans repeats
# Ou plus simple : ajouter chrZ et merger tout
echo -e "chrZ\t0\t999999999" >> repeats.bed

# Merger et trier (élimine les doublons/chevauchements automatiquement)
sort -k1,1 -k2,2n repeats.bed | bedtools merge > regions_to_exclude.bed

Ou encore plus simple, garde ton code original - bedtools subtract s'en fiche des doublons, il enlève juste tout ce qui chevauche. Donc ton code fonctionne parfaitement tel quel !

Le bedtools merge avant le subtract s'occupe déjà de fusionner les régions qui se chevauchent dans le masque de couverture de toute façon.

ok c'est pas dérangeant si des mon repeats c'est les sites qui sont pas callable car du coup mon fichier .ped c'est chr1    0    1
chr1    1    2
chr1    2    3
chr1    3    4
chr1    4    5
chr1    5    6
chr1    6    7
chr1    7    8
chr1    8    9
chr1    9    10

Ah je vois, ton fichier contient tous les sites non-callables (pas juste les répétitions), déjà au format BED !

Dans ce cas, c'est encore plus simple - ton code est parfait tel quel :
bash

# Ton fichier repeats.txt est déjà en format BED (chr start end)
# Pas besoin de conversion avec awk !

# Si tu veux quand même exclure tout chrZ en plus
echo -e "chrZ\t0\t999999999" >> repeats.txt
sort -k1,1 -k2,2n repeats.txt | bedtools merge > regions_to_exclude.bed

Ou si ton fichier repeats.txt contient déjà toutes les régions non-callables incluant chrZ, tu peux directement faire :
bash

# 3. Générer les masques .bed depuis les BAMs et enlever toutes les régions non-callables
for i in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$i]}"
    bam="${BAMS[$i]}"
    
    samtools depth -a -r ${CHR} ${bam} | \
    awk '$3 >= 5 {print $1"\t"$2-1"\t"$2}' | \
    bedtools merge | \
    bedtools subtract -a stdin -b repeats.txt > ${sample}.${CHR}.mask.bed
done

Question : Ton fichier repeats.txt contient déjà chrZ ou tu veux l'ajouter en plus ?

pour le phasage il me dit Command line: java -Xmx29694m -jar beagle.22Jul22.46e.jar
  gt=BEA_i04.chr1.vcf.gz
  out=BEA_i04.chr1.phased
  nthreads=1
No genetic map is specified: using 1 cM = 1 Mb
Reference samples:                    0
Study     samples:                    1
Window 1 [chr1:817-26529423]
Study     markers:               50,351
ERROR: there is only one sample

L'erreur vient du fait que Beagle ne peut pas phaser un seul individu sans carte génétique ou panel de référence !

Solutions possibles :

Option 1 : Utiliser WhatsHapp (meilleur pour un seul individu avec reads)
bash

module load bioinfo/whatshapp/2.2  # ou la version disponible

for i in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$i]}"
    bam="${BAMS[$i]}"
    
    whatshapp phase -o ${sample}.${CHR}.phased.vcf.gz \
        --reference genome.fasta \
        ${sample}.${CHR}.vcf.gz \
        ${bam}
done

Option 2 : Phaser tous les individus ensemble avec Beagle, puis séparer
bash

# 1. Extraire TOUS les individus ensemble
bcftools view -s CTB_1,BRE_1,BEA_1,FIR_1 -r ${CHR} ${VCF} | \
bcftools view -e 'AC==0' -O z -o all_samples.${CHR}.vcf.gz

# 2. Phaser ensemble
java -jar beagle.jar gt=all_samples.${CHR}.vcf.gz out=all_samples.${CHR}.phased

# 3. Séparer par individu après phasage
for sample in "${SAMPLES[@]}"; do
    bcftools view -s ${sample} all_samples.${CHR}.phased.vcf.gz -O z -o ${sample}.${CHR}.phased.vcf.gz
    bcftools index ${sample}.${CHR}.phased.vcf.gz
done

Option 3 : Utiliser Beagle avec carte génétique (si tu en as une)

Je recommande l'Option 2 (phaser ensemble puis séparer) - c'est le plus simple et Beagle sera plus précis avec plusieurs individus !

le problème avec l'option 2 c'est que en faisant ca je vais avoir des 0/0 quand je vais les reseparer

Ah oui exact ! Quand tu sépares après, tu vas récréer des 0/0...

Solution : phaser PUIS filtrer les 0/0
bash

# 1. Extraire TOUS les individus ensemble (GARDER les 0/0 pour l'instant)
bcftools view -s CTB_1,BRE_1,BEA_1,FIR_1 -r ${CHR} ${VCF} -O z -o all_samples.${CHR}.vcf.gz

# 2. Phaser ensemble
java -jar /usr/local/bioinfo/src/Beagle/Beagle_v5.4/beagle.22Jul22.46e.jar \
    gt=all_samples.${CHR}.vcf.gz out=all_samples.${CHR}.phased

# 3. Séparer par individu ET filtrer les 0/0 APRÈS le phasage
for sample in "${SAMPLES[@]}"; do
    bcftools view -s ${sample} all_samples.${CHR}.phased.vcf.gz | \
    bcftools view -e 'GT="0/0" || GT="0|0"' -O z -o ${sample}.${CHR}.phased.vcf.gz
    bcftools index ${sample}.${CHR}.phased.vcf.gz
done

# 4. Générer les masques .bed depuis les BAMs
for i in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$i]}"
    bam="${BAMS[$i]}"
    
    samtools depth -a -r ${CHR} ${bam} | \
    awk '$3 >= 5 {print $1"\t"$2-1"\t"$2}' | \
    bedtools merge | \
    bedtools subtract -a stdin -b repeats.txt > ${sample}.${CHR}.mask.bed
done

L'astuce : Filtrer GT="0|0" (phasé) en plus de GT="0/0" après la séparation !

Ça te va comme ça ?

ouais niquel tu peux me refaire le script entier stp

Voici le script complet :
bash

#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH -p workq
#SBATCH --mem=24G
module purge
module load bioinfo/Bcftools/1.21
module load bioinfo/Beagle/5.4
module load bioinfo/samtools/1.21
module load bioinfo/BEDTools/2.31.1

# Variables
VCF="freebayes_filt_0pct_missing_dp_filtered_masked_HWE.vcf.gz"
CHR="chr1"
REPEATS="repeats.txt"  # Ton fichier des régions non-callables
SAMPLES=("CTB_1" "BRE_1" "BEA_1" "FIR_1")
BAMS=("/work/project/loadexp/full_loadexp_data/bam/ES_Cor_i12.bam" \
      "/work/project/loadexp/full_loadexp_data/bam/FR_Bre_i01.bam" \
      "/work/project/loadexp/full_loadexp_data/bam/FR_Bea_i01.bam" \
      "/work/project/loadexp/full_loadexp_data/bam/IT_Fir_i01.bam")

# 1. Extraire TOUS les individus ensemble (garder les 0/0 pour le phasage)
bcftools view -s CTB_1,BRE_1,BEA_1,FIR_1 -r ${CHR} ${VCF} -O z -o all_samples.${CHR}.vcf.gz
bcftools index all_samples.${CHR}.vcf.gz

# 2. Phaser ensemble avec Beagle
java -jar /usr/local/bioinfo/src/Beagle/Beagle_v5.4/beagle.22Jul22.46e.jar \
    gt=all_samples.${CHR}.vcf.gz out=all_samples.${CHR}.phased

# 3. Séparer par individu ET filtrer les 0/0 APRÈS le phasage
for sample in "${SAMPLES[@]}"; do
    bcftools view -s ${sample} all_samples.${CHR}.phased.vcf.gz | \
    bcftools view -e 'GT="0/0" || GT="0|0"' -O z -o ${sample}.${CHR}.phased.vcf.gz
    bcftools index ${sample}.${CHR}.phased.vcf.gz
done

# 4. Générer les masques .bed depuis les BAMs et enlever les régions non-callables
for i in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$i]}"
    bam="${BAMS[$i]}"
    
    samtools depth -a -r ${CHR} ${bam} | \
    awk '$3 >= 5 {print $1"\t"$2-1"\t"$2}' | \
    bedtools merge | \
    bedtools subtract -a stdin -b ${REPEATS} > ${sample}.${CHR}.mask.bed
done

echo "Done! VCF phasés et masques générés pour ${CHR}"

C'est bon comme ça ?

oui merci

De rien ! Bon courage avec ton analyse MSMC2 ! 🧬

Si tu as besoin d'aide pour les étapes suivantes (génération du fichier multihetsep avec generate_multihetsep.py ou le lancement de MSMC2), n'hésite pas ! 👍

[W::vcf_parse] Contig 'chr1' is not defined in the header. (Quick workaround: index the file with tabix.)
Undefined tags in the header, cannot proceed in the sample subset mode.

Cette erreur indique que ton VCF n'a pas de ligne ##contig dans le header. Voici comment corriger :

Solution rapide :
bash

# Ajouter les contigs au header avec un fichier de référence
bcftools reheader --fai genome.fasta.fai ${VCF} -o ${VCF%.gz}_reheader.vcf.gz

# OU si tu n'as pas de .fai, créer manuellement les lignes contig
# puis utiliser la nouvelle version

Ou contourner le problème :
bash

# Au lieu de bcftools view -s pour extraire les samples
# Utilise d'abord tabix puis bcftools

# 1. Extraire le chromosome d'abord
tabix -h ${VCF} ${CHR} | bgzip > temp.${CHR}.vcf.gz
tabix -p vcf temp.${CHR}.vcf.gz

# 2. PUIS extraire les samples
bcftools view -s CTB_1,BRE_1,BEA_1,FIR_1 temp.${CHR}.vcf.gz -O z -o all_samples.${CHR}.vcf.gz

Ou la solution la plus propre - ajouter les contigs manquants :
bash

# Si tu as le .fai du génome de référence
bcftools reheader --fai reference.fasta.fai -o ${VCF%.gz}_fixed.vcf.gz ${VCF}

# Utilise ensuite le VCF fixé
VCF="${VCF%.gz}_fixed.vcf.gz"

Tu as le fichier .fai de ton génome de référence ?

nn mais c'est juste après beagle

Ah d'accord ! C'est après Beagle lors de la séparation des individus. Le problème c'est que Beagle ne copie pas toujours correctement les lignes ##contig dans le header.

Solution : Ajouter les contigs au VCF de sortie de Beagle
bash

#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH -p workq
#SBATCH --mem=24G
module purge
module load bioinfo/Bcftools/1.21
module load bioinfo/Beagle/5.4
module load bioinfo/samtools/1.21
module load bioinfo/BEDTools/2.31.1
module load bioinfo/tabix/0.2.6  # ou htslib

# Variables
VCF="freebayes_filt_0pct_missing_dp_filtered_masked_HWE.vcf.gz"
CHR="chr1"
REPEATS="repeats.txt"
SAMPLES=("CTB_1" "BRE_1" "BEA_1" "FIR_1")
BAMS=("/work/project/loadexp/full_loadexp_data/bam/ES_Cor_i12.bam" \
      "/work/project/loadexp/full_loadexp_data/bam/FR_Bre_i01.bam" \
      "/work/project/loadexp/full_loadexp_data/bam/FR_Bea_i01.bam" \
      "/work/project/loadexp/full_loadexp_data/bam/IT_Fir_i01.bam")

# 1. Extraire TOUS les individus ensemble
bcftools view -s CTB_1,BRE_1,BEA_1,FIR_1 -r ${CHR} ${VCF} -O z -o all_samples.${CHR}.vcf.gz
bcftools index all_samples.${CHR}.vcf.gz

# 2. Phaser ensemble avec Beagle
java -jar /usr/local/bioinfo/src/Beagle/Beagle_v5.4/beagle.22Jul22.46e.jar \
    gt=all_samples.${CHR}.vcf.gz out=all_samples.${CHR}.phased

# 2b. Réindexer le fichier phasé pour ajouter les contigs
tabix -p vcf all_samples.${CHR}.phased.vcf.gz

# 3. Séparer par individu ET filtrer les 0/0 APRÈS le phasage
for sample in "${SAMPLES[@]}"; do
    tabix -h all_samples.${CHR}.phased.vcf.gz ${CHR} | \
    bcftools view -s ${sample} | \
    bcftools view -e 'GT="0/0" || GT="0|0"' -O z -o ${sample}.${CHR}.phased.vcf.gz
    bcftools index ${sample}.${CHR}.phased.vcf.gz
done

# 4. Générer les masques .bed depuis les BAMs
for i in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$i]}"
    bam="${BAMS[$i]}"
    
    samtools depth -a -r ${CHR} ${bam} | \
    awk '$3 >= 5 {print $1"\t"$2-1"\t"$2}' | \
    bedtools merge | \
    bedtools subtract -a stdin -b ${REPEATS} > ${sample}.${CHR}.mask.bed
done

echo "Done! VCF phasés et masques générés pour ${CHR}"

Changement clé : Utiliser tabix -h pour extraire le chromosome avant bcftools view -s, ce qui préserve le header correctement.

réecris en anglais et met un titre au debut genre script for ...

r

# Script for ABC Random Forest parameter estimation with posterior distribution analysis

library(abcrf)
library(weights)
library(KScorrect)
library(writexl)
library(openxlsx)

# Load real data
reel <- read.table("real_data/obs_table_ref", h = T)

args <- commandArgs(trailingOnly = TRUE)
k <- as.integer(args[1])

name.list <- c("sim.Ne_ancpp", "sim.Ne_sp", "sim.Ne_wp", "sim.Ne_fu", "sim.Ne_sp_found", "sim.Ne_wp_found", "sim.Ne_fu_ancfound",
               "sim.Ne_wp_ancfound", "sim.Ne_fu_anc", "sim.Ne_wp_anc", "sim.Ne_bot_sp", "sim.Ne_bot_wp", "sim.Ne_bot_fu",
               "sim.tsplit_PP", "sim.tsplit_WP", "sim.t_bot", "sim.m_ancpp", "sim.m_spwp_anc", "sim.m_spwp_rec", "sim.m_spfu_anc",
               "sim.m_spfu_rec", "sim.m_wpfu_anc", "sim.m_wpfu_rec")

pdf(paste0("Posterior.scen_SP_WP_", name.list[k], ".pdf"), width = 300/72, height = 300/72)

predict.parameter <- data.frame("estimate_mean" = numeric(), "estimate_median" = numeric(), "90_CI_2.5" = numeric(), "90_CI_97.5" = numeric(),
                                "prior_NMAE_mean" = numeric(), "prior_NMEA_median" = numeric(),
                                "posterior_NMAE_mean" = numeric(), "posterior_NMEA_median" = numeric(),
                                "95_cov" = numeric())

df1 <- read.table("reftable_scen/reftable.recent_div.txt", h = T)

# Diagnostics: Display dimensions before selection
cat("Dimensions of real data:", dim(reel), "\n")
cat("Dimensions of reference table:", dim(df1), "\n")

reel_new <- reel[, c(1:32, 34, 36, 38:42, 44, 46, 48:52, 54, 56, 58:167)]
ref_table_1 <- as.data.frame(df1[, c(name.list[k], names(df1)[24:190])])
ref_table_1 <- ref_table_1[, c(1:33, 35, 37, 39:43, 45, 47, 49:53, 55, 57, 59:168)]
ref_table_1 <- na.omit(ref_table_1)

# Diagnostics: Check dimensions after selection
cat("Dimensions of reel_new:", dim(reel_new), "\n")
cat("Dimensions of ref_table_1:", dim(ref_table_1), "\n")
cat("Number of rows in ref_table_1 after na.omit:", nrow(ref_table_1), "\n")

# Diagnostics: Compare column names
cat("\n=== DIAGNOSTICS ===\n")
cat("Column names in ref_table_1 (excluding first column - response variable):\n")
predictors_train <- names(ref_table_1)[-1]
print(predictors_train)

cat("\nColumn names in reel_new:\n")
predictors_obs <- names(reel_new)
print(predictors_obs)

# Find missing columns
missing_in_obs <- setdiff(predictors_train, predictors_obs)
missing_in_train <- setdiff(predictors_obs, predictors_train)

if (length(missing_in_obs) > 0) {
    cat("\n!!! MISSING COLUMNS in reel_new:\n")
    print(missing_in_obs)
}

if (length(missing_in_train) > 0) {
    cat("\n!!! EXTRA COLUMNS in reel_new (not in ref_table_1):\n")
    print(missing_in_train)
}

# Correction: Ensure reel_new has exactly the same columns as ref_table_1 (except first column)
common_cols <- intersect(predictors_train, predictors_obs)
cat("\nNumber of common columns:", length(common_cols), "\n")

if (length(common_cols) < length(predictors_train)) {
    cat("\n!!! WARNING: Missing", length(predictors_train) - length(common_cols), "columns!\n")
    stop("Cannot continue: missing columns in observed data")
}

# Reorder reel_new to match ref_table_1 exactly
reel_new <- reel_new[, predictors_train]

cat("\nAfter correction - Dimensions of reel_new:", dim(reel_new), "\n")
cat("Verification: identical names?", identical(names(reel_new), predictors_train), "\n")

# Create regression model
formula <- paste("log10(", name.list[k], ") ~ .", sep = "")
formula <- as.formula(formula)
RFmodel <- regAbcrf(formula = formula,
                    data = ref_table_1,
                    ntree = 1000,
                    paral = TRUE)

# Plot 1: Prior vs posterior distribution
hist(log10(ref_table_1[, 1]),
     breaks = 100,
     col = "grey",
     ylim = c(0, 1.5),
     freq = FALSE,
     main = name.list[k],
     xlab = "",
     ylab = "probability density", cex.main = 1, cex.axis = 0.8)

posterior_RF <- predict(object = RFmodel,
                        obs = reel_new,
                        quantiles = c(0.05, 0.95),
                        training = ref_table_1,
                        paral = TRUE,
                        ncores = 32,
                        ntree = 1000,
                        rf.weights = TRUE,
                        post.err.med = TRUE)

wtd.hist(log10(ref_table_1[, 1]),
         breaks = 50,
         col = rgb(1, 1, 0, alpha = 0.5), freq = FALSE, add = TRUE,
         weight = posterior_RF$weights)

newdata <- data.frame("estimate_mean" = posterior_RF$expectation,
                      "estimate_median" = posterior_RF$med[1],
                      "90_CI_2.5" = posterior_RF$quantile[1],
                      "90_CI_97.5" = posterior_RF$quantile[2],
                      "prior_NMAE_mean" = posterior_RF$prior.NMAE.mean,
                      "prior_NMEA_median" = posterior_RF$prior.NMAE.med,
                      "posterior_NMAE_mean" = posterior_RF$post.NMAE.mean,
                      "posterior_NMEA_median" = posterior_RF$post.NMAE.med,
                      "95_cov" = posterior_RF$prior.coverage)

predict.parameter <- rbind(predict.parameter, newdata)

dev.off()

# Plot 2: Estimated vs Actual (Out-of-Bag) in PNG format
png(paste0("OOB_predictions_scen_SP_WP_", name.list[k], ".png"), 
    width = 2000, height = 2000, res = 300)

# Extract OOB predictions
oob_predictions <- RFmodel$model.rf$predictions
actual_values <- log10(ref_table_1[, 1])

# Create plot
plot(actual_values, oob_predictions,
     xlab = "Actual parameter (log10)",
     ylab = "Estimated parameter (log10)",
     main = paste("OOB predictions -", name.list[k]),
     pch = 16,
     col = rgb(0, 0, 0, alpha = 0.1),
     cex = 0.3)

# Add identity line (y=x)
abline(0, 1, col = "red", lwd = 2, lty = 2)

# Add linear regression line
abline(lm(oob_predictions ~ actual_values), col = "blue", lwd = 2)

# Calculate R²
r_squared <- cor(actual_values, oob_predictions)^2
legend("topleft", 
       legend = c(paste("R² =", round(r_squared, 3)),
                  "Identity line",
                  "Linear fit"),
       col = c("black", "red", "blue"),
       lty = c(0, 2, 1),
       lwd = c(0, 2, 2),
       pch = c(16, NA, NA),
       bty = "n",
       cex = 1.2)

dev.off()

predict.parameter$parameter <- name.list[k]
chemin_fichier_excel <- paste0("parameter_scen_SP_WP_", name.list[k], ".xlsx")
write_xlsx(predict.parameter, chemin_fichier_excel)


