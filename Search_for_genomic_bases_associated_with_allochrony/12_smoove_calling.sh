#!/bin/bash
set -euo pipefail

# Structural Variant Detection with Smoove and Duphold
# Usage: ./smoove_structural_variants.sh <sample.bam>

if [ $# -eq 0 ]; then
    echo "Usage: $0 <sample.bam>"
    exit 1
fi

# Configuration
COHORT_NAME="cohort"
REFERENCE="path/to/reference.fa"
THREADS=4
WORK_DIR=$(pwd)
SAMPLE=$1

BASE_NAME=$(basename "$SAMPLE" .bam)
BASE=$(echo "$BASE_NAME" | sed 's/\.bam$//')

mkdir -p results-smoove results-genotyped results-duphold final-vcf

# Step 1: Call SVs per sample
docker run -v ${WORK_DIR}:/work brentp/smoove \
    smoove call \
    --outdir /work/results-smoove/ \
    --name $BASE \
    --fasta /work/$REFERENCE \
    -p $THREADS \
    --genotype /work/$SAMPLE

# Step 2: Merge VCFs across samples
NUM_SAMPLES=$(ls results-smoove/*.genotyped.vcf.gz 2>/dev/null | wc -l)

if [ $NUM_SAMPLES -gt 1 ]; then
    docker run -v ${WORK_DIR}:/work brentp/smoove \
        smoove merge \
        --name merged \
        -f /work/$REFERENCE \
        --outdir /work/ \
        /work/results-smoove/*.genotyped.vcf.gz
else
    cp results-smoove/${BASE}.genotyped.vcf.gz merged.sites.vcf.gz
    tabix -p vcf merged.sites.vcf.gz
fi

# Step 3: Joint genotyping
docker run -v ${WORK_DIR}:/work brentp/smoove \
    smoove genotype \
    -d -x \
    -p $THREADS \
    --outdir /work/results-genotyped/ \
    --name ${BASE}-joint \
    --fasta /work/$REFERENCE \
    --vcf /work/merged.sites.vcf.gz \
    /work/$SAMPLE

# Step 4: Paste genotyped VCFs
NUM_GENOTYPED=$(ls results-genotyped/*.genotyped.vcf.gz 2>/dev/null | wc -l)

if [ $NUM_GENOTYPED -gt 1 ]; then
    docker run -v ${WORK_DIR}:/work brentp/smoove \
        smoove paste \
        --name ${COHORT_NAME} \
        /work/results-genotyped/*.genotyped.vcf.gz
    JOINT_VCF="${COHORT_NAME}.smoove.square.vcf.gz"
else
    JOINT_VCF="results-genotyped/${BASE}-joint.genotyped.vcf.gz"
fi

# Step 5: Annotate with duphold
docker run -v ${WORK_DIR}:/work brentp/smoove \
    duphold \
    -t $THREADS \
    -v /work/$JOINT_VCF \
    -b /work/$SAMPLE \
    -f /work/$REFERENCE \
    -o /work/results-duphold/${BASE}.duphold.vcf

bgzip -f results-duphold/${BASE}.duphold.vcf
tabix -p vcf results-duphold/${BASE}.duphold.vcf.gz

DUPHOLD_VCF="results-duphold/${BASE}.duphold.vcf.gz"

# Step 6: Filter by SV type
# Deletions: DHFFC < 0.7
bcftools view -i 'SVTYPE="DEL" && (DHFFC[0] < 0.7 || DHBFC[0] < 0.7)' \
    $DUPHOLD_VCF -Oz -o final-vcf/${BASE}.deletions.vcf.gz

# Duplications: DHBFC > 1.25
bcftools view -i 'SVTYPE="DUP" && DHBFC[0] > 1.25' \
    $DUPHOLD_VCF -Oz -o final-vcf/${BASE}.duplications.vcf.gz

# Inversions
bcftools view -i 'SVTYPE="INV"' \
    $DUPHOLD_VCF -Oz -o final-vcf/${BASE}.inversions.vcf.gz

# All filtered SVs
bcftools view -i '(SVTYPE="DEL" && DHFFC[0] < 0.7) || (SVTYPE="DUP" && DHBFC[0] > 1.25) || SVTYPE="INV"' \
    $DUPHOLD_VCF -Oz -o final-vcf/${BASE}.all_svs.filtered.vcf.gz

# Index VCFs
for vcf in final-vcf/${BASE}.*.vcf.gz; do
    tabix -p vcf $vcf
done
