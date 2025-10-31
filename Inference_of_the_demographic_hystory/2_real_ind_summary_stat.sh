#! /bin/bash
# Script for estimate individual summary statistics on real data

# VCF with individual data of LSP, LWP and FU
INPUT_VCF="demo_ind.vcf.gz"
LIST="/path/to/your/list"

mkdir -p het_ind
mkdir -p SFS
mkdir -p TajimaD

# Heterozygosity
vcftools --gzvcf $INPUT_VCF --hardy --keep ${LIST}/SP.list --out SP
awk 'NR>1' SP.hwe | awk '{split($3,el,"/"); if ((el[1]+el[2]+el[3])>0) print (el[2]/(el[1]+el[2]+el[3]))}' > het_ind/SP.het.txt
vcftools --gzvcf $INPUT_VCF --hardy --keep ${LIST}/WP.list --out WP
awk 'NR>1' WP.hwe | awk '{split($3,el,"/"); if ((el[1]+el[2]+el[3])>0) print (el[2]/(el[1]+el[2]+el[3]))}' > het_ind/WP.het.txt
vcftools --gzvcf $INPUT_VCF --hardy --keep ${LIST}/FU.list --out FU
awk 'NR>1' FU.hwe | awk '{split($3,el,"/"); if ((el[1]+el[2]+el[3])>0) print (el[2]/(el[1]+el[2]+el[3]))}' > het_ind/FU.het.txt

# SFS
vcftools --gzvcf $INPUT_VCF --counts --keep ${LIST}/SP.list --out SP
bash ind.SFS.sh SP.frq.count > SFS/SP.SFS.txt
vcftools --gzvcf $INPUT_VCF --counts --keep ${LIST}/WP.list --out WP
bash ind.SFS.sh WP.frq.count > SFS/WP.SFS.txt
vcftools --gzvcf $INPUT_VCF --counts --keep ${LIST}/FU.list --out FU
bash ind.SFS.sh FU.frq.count > SFS/FU.SFS.txt

# TajimaD
vcftools --gzvcf $INPUT_VCF --TajimaD 10000 --keep ${LIST}/SP.list --out SP
awk 'NR>1' SP.Tajima.D | awk '{print $4}' > TajimaD/SP.TajimaD.txt
vcftools --gzvcf $INPUT_VCF --TajimaD 10000 --keep ${LIST}/WP.list --out WP
awk 'NR>1' WP.Tajima.D | awk '{print $4}' > TajimaD/WP.TajimaD.txt
vcftools --gzvcf $INPUT_VCF --TajimaD 10000 --keep ${LIST}/FU.list --out FU
awk 'NR>1' FU.Tajima.D | awk '{print $4}' > TajimaD/FU.TajimaD.txt

# Create haplotype matrix
gunzip -c $INPUT_VCF | grep ^c | awk '{
  for (i = 10; i <= NF; i++) {
    split($i, field, ":");
    split(field[1], toto, "/");
    printf "\t%s\t%s", toto[1], toto[2];
  }
  print "";
}' > haplotype.ind.txt

#LSP
awk '{
  for (i = 1; i <= 50; i++) {
    printf "%s\t", $i;  
  }
  print "";  
}' haplotype.ind.txt > sp.haplo.txt

#LWP
awk '{
  for (i = 51; i <= 86; i++) {
    printf "%s\t", $i;  
  }
  print "";  
}' haplotype.ind.txt > wp.haplo.txt

#FU
awk '{
  for (i = 87; i <= 106; i++) {
    printf "%s\t", $i;  
  }
  print "";  
}' haplotype.ind.txt > fu.haplo.txt
