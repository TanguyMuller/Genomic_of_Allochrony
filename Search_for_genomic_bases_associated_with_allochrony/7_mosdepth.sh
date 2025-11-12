#! /bin/bash
# This script run mosdepth on the Z chromosome
# Command to run mosdepth : paste mosdepth/mosdepth_idsample.txt mosdepth/mosdepth_bam.txt |  awk '{print "bash mosdepth.sh", $1, $2}' | bash

id=$1
samples=$2

# Mosdepth command
mosdepth -t 8 -n --fast-mode --by 2000 ${id} ${samples}

mv ${id}.mosdepth* mosdepth/.
mv ${id}.region* mosdepth/.

# Mosdepth on the Z chromosome
gunzip -c mosdepth/${id}.mosdepth.regions.bed.gz | grep ^chrZ | gzip > Z.${id}.regions.bed.gz
