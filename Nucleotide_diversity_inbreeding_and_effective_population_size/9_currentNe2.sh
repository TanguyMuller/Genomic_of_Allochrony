#! /bin/bash
# Script for run currentNe2
#
# Command:
# while read pop; do
#  bash 9_currentNe2.sh "$pop"
# done < pop.txt

pop=$1

./currentne2 -r 4 -x -t 16 ${pop}.pruned_mono.ped -o ${pop}_x_all
