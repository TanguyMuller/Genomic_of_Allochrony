### Script for variant calling on chrZ

#! /bin/bash

# Set input variables
BAMLIST=$1
REGION=$2
PATH_TO_ASSEMBLY=/path/to/assembly/PPM.fa
base_name=$(basename "$BAMLIST" .portugal.bam.list)

# FreeBayes is slower than other variant callers
# It is run on 1 Mb non-overlapping regions

# Variant calling for Z chromosome (haploid samples defined in sample.bed)
# Note: -A assigns ploidy 1 to listed samples (e.g., females on chrZ)
freebayes -f $PATH_TO_ASSEMBLY -L $BAMLIST -r $REGION \
  -E -1 --max-complex-gap -1 --haplotype-length -1 -K \
  -C 1 -F 0.01 -G 1 --limit-coverage 500 -n 3 -m 30 -q 20 \
  -A sample.bed | \
gzip -c > res_vcf/${base_name}.${REGION}.chrZ.freebayes.vcf.gz

### Command to run 4_variant_calling_with_freebayes_chrZ.sh

# Note : the format of the regions.chunks file is chr:start-end
#awk '{print "sbatch -J FBch"FNR" --output FBch"FNR" 4_variant_calling_with_freebayes_chrZ.sh all.portugal.bam.list",$1}' regions.Z.chunks
#awk '{print "sbatch -J FBch"FNR" --output FBch"FNR" 4_variant_calling_with_freebayes_chrZ.sh pool.portugal.bam.list",$1}' regions.Z.chunks
#awk '{print "sbatch -J FBch"FNR" --output FBch"FNR" 4_variant_calling_with_freebayes_chrZ.sh ind.portugal.bam.list",$1}' regions.Z.chunks
