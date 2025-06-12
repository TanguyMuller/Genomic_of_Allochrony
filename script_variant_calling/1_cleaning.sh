### Script for cleaning raw fastq

#!/bin/bash

# Set input variables
QUERY_NAME=$1
NTHREADS=4

# Cleaning raw fastq files
fastp -i $QUERY_NAME.fq_1.gz -I $QUERY_NAME.fq_2.gz -o $QUERY_NAME.clean.fq_1.gz -O $QUERY_NAME.clean.fq_2.gz \
-w $NTHREADS -h $QUERY_NAME.html -j $QUERY_NAME.jston

rm $QUERY_NAME.fq_1.gz $QUERY_NAME.fq_2.gz

### Command to run cleaning.sh

#ls *fq_1.gz |\
#awk '{split($1,tmp,"\.");{print "sbatch -J",tmp[1],"--output",tmp[1],"1_cleaning.sh",tmp[1]}}' - |\
#bash
