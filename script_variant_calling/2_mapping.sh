### Script for mapping, filtering q20, sorting and bringing stats

#!/bin/bash

# Set input variables
PATH_TO_ASSEMBLY=/path/to/assembly/PPM.fa
QUERY_NAME=$1
STUDY_NAME=$2
NTHREADS=4
RGLB='LIB-'$STUDY_NAME
RGSM=$STUDY_NAME
SAMTOOLS_THREADS=$(($NTHREADS-1))

# Mapping reads to the reference genome with BWA-MEM2
# Adds read group information and filters alignments with mapping quality < 20
bwa-mem2 mem -M -t $NTHREADS -R @RG\\tID:${STUDY_NAME}\\tPL:ILLUMINA\\tLB:${RGLB}\\tSM:${RGSM} \
$PATH_TO_ASSEMBLY $QUERY_NAME.clean.fq_1.gz $QUERY_NAME.clean.fq_2.gz | \
samtools view --threads $SAMTOOLS_THREADS -q 20 -uS - | \

# Sorting alignments by read name (required for fixmate)
samtools sort -T $STUDY_NAME -n --threads $SAMTOOLS_THREADS -O BAM -o $STUDY_NAME.q20.sorted.bam -

# Fixing mate information
samtools fixmate -rpcm --threads $SAMTOOLS_THREADS $STUDY_NAME.q20.sorted.bam $STUDY_NAME.q20.sorted.fm.bam
rm $STUDY_NAME.q20.sorted.bam

# Sorting alignments by coordinate (required before markdup)
samtools sort -T $STUDY_NAME --threads $SAMTOOLS_THREADS $STUDY_NAME.q20.sorted.fm.bam > $STUDY_NAME.q20.sorted.bam
rm $STUDY_NAME.q20.sorted.fm.bam

# Computing alignment statistics (pre-deduplication)
samtools stats --threads $SAMTOOLS_THREADS $STUDY_NAME.q20.sorted.bam > $STUDY_NAME.q20.bamstats

# Marking and removing duplicate reads
samtools markdup --threads $SAMTOOLS_THREADS -r $STUDY_NAME.q20.sorted.bam $STUDY_NAME.q20.dedup.bam

# Indexing the deduplicated BAM file
samtools index -@ $SAMTOOLS_THREADS $STUDY_NAME.q20.dedup.bam

# Computing alignment statistics (post-deduplication)
samtools stats --threads $SAMTOOLS_THREADS $STUDY_NAME.q20.dedup.bam > $STUDY_NAME.q20.dedup.bamstats

rm $STUDY_NAME.q20.sorted.bam

### Command to run 2_mapping.sh

#ls *fq_1.gz |\
#awk 'NR==FNR { names[NR]=$0; next } {split($1,tmp,"."); print "sbatch -J",tmp[1],"--output",tmp[1],"2_mapping.sh",tmp[1],names[FNR] }' sample_portugal.list - |\
#bash
