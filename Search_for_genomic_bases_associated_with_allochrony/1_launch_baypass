#!/bin/bash

baypass="/path/to/baypass/baypass"
num=$1
outpref='AsubAcRc_'${num}

# the autosomes baypass analyses was divided in 100 sub sample
$baypass -poolsizefile taille_haploide_pool_auto.txt -pooldatafile A.sub.rcfile.${num} -countdatafile A.sub.acfile.${num} -contrastfile contrast -outprefix $outpref
bzip2 ${outpref}*out
gzip A.sub.?cfile.${num}

$baypass -poolsizefile taille_haploide_pool_Z.txt -pooldatafile Z.baypass.rcfile -gldatafile Z.baypass.glfile -contrastfile contrast -outprefix Z_gl -nthreads 8
bzip2 Z_gl*out
