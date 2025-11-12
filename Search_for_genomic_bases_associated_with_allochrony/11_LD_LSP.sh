#!/bin/bash
# Script for have the LD by block for the highly differenciated regions on the Z chromosome for LSP population

LDBlockShow -InVCF LSP.vcf.gz -SeleVar 2 -OutPut LSP_13500000-15000000 -Region chrZ:13500000-15000000
