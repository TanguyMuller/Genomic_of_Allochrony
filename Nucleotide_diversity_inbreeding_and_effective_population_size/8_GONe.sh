#!/bin/bash
# Script to run GONE on 100 independent subsets of 50,000 randomly
# selected SNPs. This version is designed for use on a cluster
# (e.g., SLURM-based) to parallelize runs efficiently.
#
# Commands:
# for j in {1..100}; do awk -v i="$j" '{print "sbatch 8_GONe.sh", $1, "hc001", i}' pop.txt; done > command_freebayes_hc001
# for j in {1..100}; do awk -v i="$j" '{print "sbatch 8_GONe.sh", $1, "hc005", i}' pop.txt; done > command_freebayes_hc005
#
# bash command_freebayes_hc001
# bash command_freebayes_hc005

pop=$1
hc=$2
analyse=$3
mkdir -p ${hc}/${pop}
path_vcf="/path/to/pop/vcf/map/and/ped/"

# PROGRAMMES, script_GONE.sh and INPUT_PARAMETERS_FILE can be found in the github of GONe: https://github.com/esrud/GONE
# Note: the parameter hc was put at 0.01 or 0.05 in the INPUT_PARAMETERS_FILE
cp -r PROGRAMMES script_GONE.sh INPUT_PARAMETERS_FILE ${hc}/${pop}/.

# Create analysis directory
cd ${hc}/${pop}
mkdir analyses_${analyse}
cd analyses_${analyse}

# Run GONE via SLURM
# Each subset is submitted as a separate SLURM job.
cp -r ../PROGRAMMES ../script_GONE.sh ../INPUT_PARAMETERS_FILE .
sbatch script_GONE.sh ${path_vcf}/${pop}.pruned_mono
