#!/bin/bash
# Create an mpileup file of the three outgroups restricted to the positions present in the VCF.
# Command to run : bash 1_create_mpileup.sh

PATH_TO_ASSEMBLY="/your/path/to/fasta"
PATH_TO_FOLDER="/your/path/to/folder"

mkdir -p "${PATH_TO_FOLDER}/mpileup"
mkdir -p "${PATH_TO_FOLDER}/list"

# Extract positions (CHROM, POS) from the VCF
gunzip -c "${PATH_TO_FOLDER}/vcf/demo.portugal.vcf.gz" | grep -v "^#" | awk '{print $1"\t"$2}' > "${PATH_TO_FOLDER}/list/positions_load_vcf.txt"

# Create mpileup for outgroups restricted to VCF positions
samtools mpileup -B --ignore-RG -f "${PATH_TO_ASSEMBLY}" -l "${PATH_TO_FOLDER}/list/positions_load_vcf.txt" -b list_ind_bam_outgroup.txt > "${PATH_TO_FOLDER}/mpileup/load_outgroup.mpileup"
