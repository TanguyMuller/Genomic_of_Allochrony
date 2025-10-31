#!/bin/bash
# Script for run GONe 100 independent subsets of 50,000 randomly selected SNPs
# We run this in genobioinfo cluster for parallelize

VCF=$1
pop=$2
hc=$3
analyse=$4
mkdir -p ${hc}/${pop}
path_vcf="/path/to/pop/vcf/map/and/ped/"

# PROGRAMMES, script_GONE.sh and INPUT_PARAMETERS_FILE can be found in the github of GONe: https://github.com/esrud/GONE
# Note: the parameter hc was put at 0.01 or 0.05 in the INPUT_PARAMETERS_FILE
cp -r PROGRAMMES script_GONE.sh INPUT_PARAMETERS_FILE ${hc}/${pop}/.

cd ${hc}/${pop}
mkdir analyses_${analyse}
cd analyses_${analyse}

cp -r ../PROGRAMMES ../script_GONE.sh ../INPUT_PARAMETERS_FILE .
sbatch script_GONE.sh ${path_vcf}/${pop}_10covperpool
