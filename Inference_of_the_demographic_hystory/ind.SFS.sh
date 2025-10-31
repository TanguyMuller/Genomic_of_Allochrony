#!/bin/bash
# Script for estimate SFS for each population from simulate individuals data.
# This script is runned by ind.sh

# Initialize an associative array to hold counts of each allele number
declare -A counts

input_file=$1

# Get the max allele count 
max_alleles=$(awk 'NR > 1 {print $4; exit}' "$input_file")
awk -v max_alleles="$max_alleles" -F '\t' '
BEGIN {
    # Initialize counters for each possible allele count (0 to max_alleles)
    for (i = 0; i <= max_alleles; i++) {
        count[i] = 0
    }
}
NR > 1 { 
    allele_counts = $6         
    gsub(/[{}]/, "", allele_counts)  
    split(allele_counts, ac, ",")    
    
    for (i in ac) {
        split(ac[i], pair, ":")      
        allele_count = pair[2]       
        count[allele_count]++        
    }
}
END {
    for (i = 0; i <= max_alleles; i++) {
        print count[i]
    }
}' "$input_file"
