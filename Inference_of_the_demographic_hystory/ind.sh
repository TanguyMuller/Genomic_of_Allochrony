#! /bin/bash
# Script runned by concomitent_div.R, recent_div.R or old_div.R for estiamate summary statistics from simulate individual data 

seed=$1

# Concatenate 1000 vcf create from simu.py for each simulation
ls all_chromosomes_seed_$seed* > list_vcf_$seed
bcftools concat -f list_vcf_$seed -Oz > $seed.vcf.gz
rm all_chromosomes_seed_$seed*
rm list_vcf_$seed

# For 25 LSP individuals 
for i in {0..24};do
echo "tsk_"$i >> list_SP_$seed
done

# For 18 LWP individuals
for i in {25..42};do
echo "tsk_"$i >> list_WP_$seed
done

# For 10 FU individuals
for i in {43..52};do
echo "tsk_"$i >> list_FU_$seed
done

# Create population vcf
bcftools view -S list_SP_$seed $seed.vcf.gz --force-samples -Oz > SP_$seed.vcf.gz
bcftools view -S list_WP_$seed $seed.vcf.gz --force-samples -Oz > WP_$seed.vcf.gz
bcftools view -S list_FU_$seed $seed.vcf.gz --force-samples -Oz > FU_$seed.vcf.gz

# Heterozygosity
vcftools --gzvcf SP_$seed.vcf.gz --hardy --out SP_$seed
awk 'NR>1' SP_$seed.hwe | awk '{split($3,el,"/"); if ((el[1]+el[2]+el[3])>0) print (el[2]/(el[1]+el[2]+el[3]))}' > SP_$seed.txt
vcftools --gzvcf WP_$seed.vcf.gz --hardy --out WP_$seed
awk 'NR>1' WP_$seed.hwe | awk '{split($3,el,"/"); if ((el[1]+el[2]+el[3])>0) print (el[2]/(el[1]+el[2]+el[3]))}' > WP_$seed.txt
vcftools --gzvcf FU_$seed.vcf.gz --hardy --out FU_$seed
awk 'NR>1' FU_$seed.hwe | awk '{split($3,el,"/"); if ((el[1]+el[2]+el[3])>0) print (el[2]/(el[1]+el[2]+el[3]))}' > FU_$seed.txt

# SFS
vcftools --gzvcf SP_$seed.vcf.gz --counts --out SP_$seed
bash ind.SFS.sh SP_$seed.frq.count > SP_$seed.SFS.txt
vcftools --gzvcf WP_$seed.vcf.gz --counts --out WP_$seed
bash ind.SFS.sh WP_$seed.frq.count > WP_$seed.SFS.txt
vcftools --gzvcf FU_$seed.vcf.gz --counts --out FU_$seed
bash ind.SFS.sh FU_$seed.frq.count > FU_$seed.SFS.txt

# Tajima D
vcftools --gzvcf SP_$seed.vcf.gz --TajimaD 10000 --out SP_$seed
awk 'NR>1' SP_$seed.Tajima.D | awk '{print $4}' > SP_$seed.TajimaD.txt
vcftools --gzvcf WP_$seed.vcf.gz --TajimaD 10000 --out WP_$seed
awk 'NR>1' WP_$seed.Tajima.D | awk '{print $4}' > WP_$seed.TajimaD.txt
vcftools --gzvcf FU_$seed.vcf.gz --TajimaD 10000 --out FU_$seed
awk 'NR>1' FU_$seed.Tajima.D | awk '{print $4}' > FU_$seed.TajimaD.txt

chmod +x *.txt

rm *_$seed.hwe
rm *_$seed.frq.count
rm *_$seed.log
rm *_$seed.Tajima.D
rm list_*_$seed
rm $seed.vcf.gz
rm *_$seed.vcf.gz
