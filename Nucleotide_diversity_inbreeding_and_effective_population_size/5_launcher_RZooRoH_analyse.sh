#!/bin/bash
# Script for run RZooRoH analyses for each population
# Note: launch this in parallel in slurm environemen
#
# Command
# bash r_launcher_RZooRoH_analyse.sh 

dataset=$1

echo "RZooRoH analyses on $dataset dataset"
awk -v vcf="$dataset" '{print "sbatch RZooRoH_analyse.sh",vcf,$1}' pop.txt | bash
echo "RZooRoH is running for $dataset"

### Command
