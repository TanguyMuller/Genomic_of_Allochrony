### Script for filt raw vcf files 

#!/bin/bash -l

# setting vcf file name
vcf_file=$1
chr=$2
base_name=$(basename "$vcf_file" .vcf)
    
# Defining the output file name
output_file="${base_name}_filt"

# Initializing statistics file
echo "FILTERING STATISTICS FOR ${base_name}" >> count.snp.${chr}.txt
echo "----------------------------------------" >> count.snp.${chr}.txt
echo "Date: $(date)" >> count.snp.${chr}.txt
echo "" >> count.snp.${chr}.txt
    
# Running the bcftools commands
echo "Decomposition and initial filtering..." >> count.snp.${chr}.txt
vt decompose -s "$vcf_file" -o tmp.${base_name}.vcf.gz 
bcftools index tmp.${base_name}.vcf.gz
vt decompose_blocksub -a tmp.${base_name}.vcf.gz | \
bcftools filter -i 'QUAL>=30 & TYPE="SNP" & INFO/AO >= 4 & INFO/AO <= (INFO/DP-4) & SAF > 0 & SAR > 0 & RPR > 0 & RPL > 0 & (STRLEN(REF)==1) & (STRLEN(ALT)==1)' --SnpGap 2 -Ov - | \
bcftools norm -m + -Ov - | bcftools view -m2 -M2 -V indels -Oz > tmp.$output_file.vcf.gz
bcftools index tmp.$output_file.vcf.gz
rm tmp.${base_name}.vcf.gz*

echo "Number of bialleliq SNP with quality >=30 : $(bcftools view -H tmp.$output_file.vcf.gz | wc -l)" >> count.snp.${chr}.txt

# The decomposition command introduces hemizygous genotypes
bcftools filter -e 'GT="0/." || GT="./1" || GT="./."' -Oz tmp.$output_file.vcf.gz > tmp.${output_file}_half.vcf.gz
bcftools index tmp.${output_file}_half.vcf.gz
echo "Number of SNPs after excluding genotypes with missing alleles (GT=0/.|GT=./1|GT=./., but keeping fully missing noted GT=.): $(bcftools view -H tmp.${output_file}_half.vcf.gz | wc -l)" >> count.snp.${chr}.txt
rm tmp.$output_file.vcf.gz*

# Reducing VCF file size by removing unneeded INFO fields: this step is not really necessary and takes a lot of computing time
echo "" >> count.snp.${chr}.txt
echo "Reducing VCF file size by removing unneeded INFO fields..." >> count.snp.${chr}.txt
gunzip -c tmp.${output_file}_half.vcf.gz | \
gawk ' 
BEGIN { 
    OFS = "\t"; 
    ORS = ""; 
    ENVIRON["LC_ALL"] = "C"; 
} 
{ 
    # Handling VCF headers
    if (substr($1, 1, 1) == "#") { 
        if (/fileformat/ || /contig/ || /ID=RO/ || /ID=AO/ || /ID=GT,/ || /ID=GQ/ || /ID=GL/ || /ID=DP,/ || /ID=AD/) { 
            print $0 "\n"; 
        } 
        if (/CHROM/) { 
            print $0 "\n"; 
            header_found = 1; 
        } 
    } 
    # Individual data
    else if (header_found) { 
        # Initialisation des indices de champs
        if (!fields_initialized) { 
            nf = split($9, tmp, ":"); 
            # Initialiser tous les champs à -1 (non trouvé)
            gt_field = dp_field = ad_field = gl_field = ro_field = ao_field = -1;
            
            for (i = 1; i <= nf; i++) { 
                if (tmp[i] == "GT") gt_field = i; 
                if (tmp[i] == "DP") dp_field = i; 
                if (tmp[i] == "AD") ad_field = i; 
                if (tmp[i] == "GL") gl_field = i; 
                if (tmp[i] == "RO") ro_field = i; 
                if (tmp[i] == "AO") ao_field = i; 
            } 
            out_fmt_field = "GT:DP:AD:RO:AO:GL"; 
            fields_initialized = 1; 
        } 
        
        # Sum of DP for all individuals
        dp_sum = 0; 
        for (j = 10; j <= NF; j++) { 
            split($j, tmp, ":"); 
            
            if (dp_field > 0 && length(tmp) >= dp_field && tmp[dp_field] != "." && tmp[dp_field] != "") {
                dp_sum += tmp[dp_field]; 
            }
        } 
        
        # Extracting and reorganizing the INFO field
        split($8, info, ";");
        
        info_str = "";
        if (length(info) >= 6) {
            info_str = info[6];
        }
        
        info_ab = "";
        if (length(info) >= 29) {
            info_ab = info[29];
        }
        
        info_out = info_str;
        if (info_str != "") {
            info_out = info_out ";";
        }
        info_out = info_out "DP=" dp_sum;
        if (info_ab != "") {
            info_out = info_out ";" info_ab;
        }
        
        print $1, $2, $3, $4, $5, $6, ".", info_out, out_fmt_field; 
        
        # Extraction of individual data
        for (j = 10; j <= NF; j++) { 
            split($j, sample_data, ":"); 
            
            # Préparation des valeurs de sortie avec gestion des cas manquants
            gt_val = ".";
            dp_val = ".";
            ad_val = ".";
            ro_val = ".";
            ao_val = ".";
            gl_val = ".";
            
            if (gt_field > 0 && length(sample_data) >= gt_field && sample_data[gt_field] != "") {
                gt_val = sample_data[gt_field];
            }
            
            if (dp_field > 0 && length(sample_data) >= dp_field && sample_data[dp_field] != "") {
                dp_val = sample_data[dp_field];
            }
            
            if (ad_field > 0 && length(sample_data) >= ad_field && sample_data[ad_field] != "") {
                ad_val = sample_data[ad_field];
            }
            
            if (ro_field > 0 && length(sample_data) >= ro_field && sample_data[ro_field] != "") {
                ro_val = sample_data[ro_field];
            }
            
            if (ao_field > 0 && length(sample_data) >= ao_field && sample_data[ao_field] != "") {
                ao_val = sample_data[ao_field];
            }
            
            if (gl_field > 0 && length(sample_data) >= gl_field && sample_data[gl_field] != "") {
                gl_val = sample_data[gl_field];
            }
            
            print "\t" gt_val ":" dp_val ":" ad_val ":" ro_val ":" ao_val ":" gl_val;
        } 
        print "\n"; 
    } 
}' | bcftools view -Oz > tmp.${output_file}_half_reduce.vcf.gz

# Using bcftools for index the file
bcftools index tmp.${output_file}_half_reduce.vcf.gz

echo "" >> count.snp.${chr}.txt
echo "FILTERING ON DEPTH " >> count.snp.${chr}.txt
echo "----------------------------------------------------" >> count.snp.${chr}.txt
echo "" >> count.snp.${chr}.txt

# Number of individuals in the VCF file
total_individuals_pools=$(bcftools query -l tmp.${output_file}_half_reduce.vcf.gz | wc -l)
echo "Number of individuals: $total_individuals_pools" >> count.snp.${chr}.txt

# Fitering on extremely cover sites and missing genotypes
echo "" >> count.snp.${chr}.txt

if [[ "$base_name" == "all.portugal.freebayes" | "$base_name" == "all.portugal.freebayes.chrZ" ]]; then
	echo "Filtering on extremely cover sites and missing genotypes..." >> count.snp.${chr}.txt
	bcftools view -i "F_PASS(FORMAT/DP>0) == 1 & INFO/DP>200 & INFO/DP<5000" tmp.${output_file}_half_reduce.vcf.gz -Oz -o ${output_file}.vcf.gz
	bcftools index ${output_file}_0pct_missing_dp_filtered.vcf.gz
	echo "Number of SNP after filtering on depth : $(bcftools view -H ${output_file}_0pct_missing_dp_filtered.vcf.gz | wc -l)" >> count.snp.${chr}.txt
else
	echo "Filtering on extremely cover sites and missing genotypes..." >> count.snp.${chr}.txt
	bcftools view -i "F_PASS(FORMAT/DP>0) == 1 & INFO/DP>200 & INFO/DP<4000" tmp.${output_file}_half_reduce.vcf.gz -Oz -o ${output_file}.vcf.gz
	bcftools index ${output_file}_0pct_missing_dp_filtered.vcf.gz
	echo "Number of SNP after filtering on depth : $(bcftools view -H ${output_file}_0pct_missing_dp_filtered.vcf.gz | wc -l)" >> count.snp.${chr}.txt
fi

# Removing temporary files
rm -f tmp.${output_file}_half.vcf.gz* tmp.${output_file}_half_reduce.vcf.gz* 

# Rename the file for simplicity
mv ${output_file}_0pct_missing_dp_filtered.vcf.gz ${output_file}.vcf.gz

echo "End filtering"

### Command to run 6_filtering_vcf.sh

#bash 6_filtering_vcf.sh all.portugal.freebayes.vcf.gz Auto
#bash 6_filtering_vcf.sh pool.portugal.freebayes.vcf.gz Auto
#bash 6_filtering_vcf.sh ind.portugal.freebayes.vcf.gz Auto
#bash 6_filtering_vcf.sh all.portugal.freebayes.chrZ.vcf.gz chrZ
#bash 6_filtering_vcf.sh all.portugal.freebayes.chrZ.vcf.gz chrZ
#bash 6_filtering_vcf.sh all.portugal.freebayes.chrZ.vcf.gz chrZ
