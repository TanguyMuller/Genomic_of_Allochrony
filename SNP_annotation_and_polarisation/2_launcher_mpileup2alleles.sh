#!/bin/bash
# Create a table of alleles from an mpileup file for the three outgroups.
# Command to run : bash 2_launcher_mpileup2alleles.sh

perl script/mpileup2allelesOutgroup_extract_3_og.pl \
  mpileup/load_outgroup.mpileup \
  mpileup/load_outgroup.txt

perl script/mpileup2allelesOutgroup_extract_3_og.pl \
  mpileup/Z_load_outgroup.mpileup \
  mpileup/Z_load_outgroup.txt
