#! /bin/bash
# Script for create the reference table
# Script on slurm

#SBATCH --mem=32G
#SBATCH --export=ALL
#SBATCH -c 1
#SBATCH -n 1


##SBATCH --mem=64G #memory reservation
##SBATCH --time=00:10:00 #job time limit
##SBATCH --mail-user=[address] 
##SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH -J testjob #job name
##SBATCH -o output.out #output file name
##SBATCH -e error.out #error file name


## Re assignate VARS

function error {
	echo -e 1>&2 "Usage: $0 seed nb_inc\n\t$1"
	exit 1
}

[ $# == 2 ] || error "bad parameters"
[[ "$1" =~ ^[[:digit:]]+$ ]] || error "seed must be integer, bad value $1"
[[ "$2" =~ ^[[:digit:]]+$ ]] || error "nb_inc must be integer, bad value $1"

seed_env=$1
steps=$2

module purge
module load devel/python/Python-3.7.9
module load statistics/R/4.3.0

echo "Simul with seed between $seed_env and $((seed_env + steps -1)), "
mkdir dr_$seed_env
cp *.sh *.R *.py dr_$seed_env/.
cd dr_$seed_env/.
Rscript *.R $seed_env $steps
echo "---thats all folks"
