#!/bin/bash -l
#SBATCH --ntasks 1
#SBATCH --mem=32G
# Script to run RZooROH analyses. Need to run 5_launcher_RZooRoH_analyse.sh to run this code.
# Note: This code was run on genobioinfo cluster 


# set sofware genobioinfo
module purge
module load statistics/R/4.3.0
module load bioinfo/VCFtools/0.1.16
module load bioinfo/Bcftools/1.17
module load bioinfo/PLINK/1.90b7

dataset=$1
pop=$2
path_list="/your/path/to/list"
path_script="/your/path/to/script"
base_name=$(basename "${dataset}" .vcf.gz)
mkdir -p vcf

# Create a VCF for each population
bcftools view -S ${path_list}/${pop}_ind.list -T ${path_list}/positions_${base_name}_10covperpool.txt --force-samples -Oz ${dataset} > vcf/${pop}_${base_name}.vcf.gz

# Convert VCF file in oxford format for RZooRoH
plink --vcf vcf/${pop}_${base_name}.vcf.gz --recode --out vcf/${pop}_${base_name} --double-id --chr-set 49 --allow-extra-chr
plink --file vcf/${pop}_${base_name} --recode oxford --out vcf/${pop}_${base_name} --double-id --chr-set 49 --allow-extra-chr

# GL for each individual for each population to convert GL in phred score in R for RZooRoh
gunzip -c vcf/${pop}_${base_name}.vcf.gz | grep ^c | awk 'BEGIN {OFS="\t"}
{
    printf "%s\t%s", $1, $2;
    for (i=10; i<=NF; i++) {
        split($i, tmp, ":"); 
        split(tmp[6], toto, ",");
        printf "\t%s\t%s\t%s", toto[1], toto[2], toto[3];
    }
    print "";  # This prints a newline at the end of each record
}' > GL/${pop}_${base_name}_GL.txt


# Convert GL data in PL data
echo "Conversion begins for ${pop}"
Rscript ${path_script}/convert_gl_to_pl.R ${pop}_${base_name}
echo "Conversion end for ${pop}"

# RZooRoH analyses
echo "RZooRoH begins for ${pop}"
Rscript ${path_script}/RZooRoH.R $pop ${base_name}
echo "RZooRoH end for ${pop}"
