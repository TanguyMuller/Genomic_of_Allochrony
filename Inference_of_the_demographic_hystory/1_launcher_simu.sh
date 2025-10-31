#! /bin/bash
# Script for create the reference table
# This version is designed for use on a cluster
# (e.g., SLURM-based) to parallelize runs efficiently.
#
# Command to run 50,000 simulations in 100 batches of 500 each:
# for i in $(seq 1 500 49501); do
#   sbatch create_ref_table.sh $i 500
# done

#SBATCH --mem=32G
#SBATCH --export=ALL
#SBATCH -c 1
#SBATCH -n 1

seed_env=$1
steps=$2

# Load required modules
module purge
module load devel/python/Python-3.7.9
module load statistics/R/4.3.0

# Create a dedicated directory for this simulation batch
echo "Running simulations with seeds between $seed_env and $((seed_env + steps - 1))"
mkdir dr_$seed_env

# Copy necessary scripts to the working directory
cp *.sh *.R *.py dr_$seed_env/
cd dr_$seed_env/

# Run the R simulation script with the provided seed range
Rscript *.R $seed_env $steps

echo "--- All simulations completed successfully ---"

