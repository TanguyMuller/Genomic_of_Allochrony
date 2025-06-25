### Script for variant calling 

#! /bin/bash

# Set input variables
BAMLIST=$1
REGION=$2
PATH_TO_ASSEMBLY=/path/to/assembly/PPM.fa
base_name=$(basename "$BAMLIST" .portugal.bam.list)

# FreeBayes is slower than other variant callers
# It is run on 1 Mb non-overlapping regions

# Variant calling for autosomal regions
freebayes -f $PATH_TO_ASSEMBLY -L $BAMLIST -r $REGION \
  -E -1 --max-complex-gap -1 --haplotype-length -1 -K \
  -C 1 -F 0.01 -G 1 --limit-coverage 500 -n 3 -m 30 -q 20 | \
gzip -c > res_vcf/${base_name}.${REGION}.freebayes.vcf.gz

### Command to run 3_variant_calling_with_freebayes.sh

# Note : the format of the regions.chunks file is chr:start-end
#awk '{print "sbatch -J FBch"FNR" --output FBch"FNR" 3_variant_calling_with_freebayes.sh all.portugal.bam.list",$1}' regions.chunks
#awk '{print "sbatch -J FBch"FNR" --output FBch"FNR" 3_variant_calling_with_freebayes.sh pool.portugal.bam.list",$1}' regions.chunks
#awk '{print "sbatch -J FBch"FNR" --output FBch"FNR" 3_variant_calling_with_freebayes.sh ind.portugal.bam.list",$1}' regions.chunks
