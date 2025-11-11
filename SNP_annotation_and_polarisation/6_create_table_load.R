# Script to prepare input data for genetic load analysis.

library(tidyverse)
library(dplyr)
library(openxlsx)
library(data.table)

# Load VCF table (genotypes + SnpEff annotations)
tab <- read.delim("vcf/load.ann.GT.txt", skip = 1, header = FALSE)

# Define population labels
SP <- paste0("SP_", 1:25)
WP <- paste0("WP_", 1:18)
FU <- paste0("FU_", 1:10)
pop <- c(SP, WP, FU)
det.pop <- c("SP", "WP", "FU")

# Rename columns
colnames(tab) <- c("Chromosome", "position", "REF", "ALT", SP, WP, FU, "SNPEFF_2", "SNPEFF_3", "SNPEFF_11")

# Convert relevant columns to factors
tab$SNPEFF_2 <- as.factor(tab$SNPEFF_2)
tab$SNPEFF_3 <- as.factor(tab$SNPEFF_3)
tab$SNPEFF_11 <- as.factor(tab$SNPEFF_11)

# Add SNP ID and annotation columns
tab$SNPbi <- paste(tab$Chromosome, "__", tab$position, sep = "")
tab$AA1 <- as.factor(str_sub(tab$SNPEFF_11, 3, 5))
tab$AA2 <- as.factor(ifelse(tab$SNPEFF_3 == "HIGH", NA, str_sub(tab$SNPEFF_11, -3, -1)))
tab$exon_position_status <- as.factor(ifelse(is.na(tab$SNPEFF_11), "not_exon", "EXON"))
tab$intron_position_status <- as.factor(ifelse(tab$SNPEFF_2 == "intron_variant", "INTRON", "not_intron"))
tab$intergenic_position_status <- as.factor(ifelse(tab$SNPEFF_2 == "intergenic_region", "INTERGENIC", "not_intergenic"))

cat(paste0("Initial number of SNPs: ", nrow(tab), "\n"), file = "statistiques_load.txt", append = TRUE)

################################################
# Load and merge outgroup data
################################################

outgroup <- read.delim("polarized_result.txt")[, c(1, 14, 15)]
tab_all <- merge(tab, outgroup, by = "SNPbi", all = FALSE)
rm(tab); gc()

################################################
# Polarization and allele correction
################################################

tab_pola <- tab_all %>%
  mutate(across(6:58, as.character)) %>%
  mutate(across(6:58, ~ case_when(
    REF_is_ANC == 1 & . == "0/0" ~ "0/0",
    REF_is_ANC == 0 & . == "0/0" ~ "1/1",
    . == "0/1" ~ "0/1",
    REF_is_ANC == 0 & . == "1/1" ~ "0/0",
    REF_is_ANC == 1 & . == "1/1" ~ "1/1",
    TRUE ~ .
  )))

rm(tab_all); gc()

# Keep only sites where ancestral allele matches REF or ALT
tab_pola <- tab_pola[tab_pola$ancestral_allele %in% c(tab_pola$REF, tab_pola$ALT), ]

# Swap REF and ALT when REF is not ancestral
swap_indices <- which(tab_pola$REF_is_ANC == 0)
tab_pola$ALT[swap_indices] <- tab_pola$REF[swap_indices]
tab_pola$REF[swap_indices] <- tab_pola$ancestral_allele[swap_indices]
tab_pola <- tab_pola[, -c(174, 175)]

################################################
# Split genotype fields (0/0, 0/1, 1/1)
################################################

setDT(tab_pola)
split_genotype <- function(geno_vec) {
  alleles <- tstrsplit(geno_vec, "/", type.convert = TRUE)
  alleles[[1]][is.na(geno_vec)] <- NA
  alleles[[2]][is.na(geno_vec)] <- NA
  return(alleles)
}

for (col_name in pop) {
  cat("Processing", col_name, "\n")
  geno <- tab_pola[[col_name]]
  split <- split_genotype(geno)
  tab_pola[[paste0(col_name, "_allele_1")]] <- split[[1]]
  tab_pola[[paste0(col_name, "_allele_2")]] <- split[[2]]
}

# Filter invalid amino acids
tab_pola <- subset(tab_pola, !(AA1 == "???") | is.na(AA1))
tab_pola <- subset(tab_pola, !(AA2 == "t1?") | is.na(AA2))

cat(paste0("Number of well-polarized SNPs: ", nrow(tab_pola), "\n"),
    file = "statistiques_load.txt", append = TRUE)

################################################
# Compute allele counts and genotype categories
################################################

tab_pola <- as.data.frame(tab_pola)

for (p in det.pop) {
  cols <- grep(paste0("^", p, "_[0-9]+_allele_[12]$"), names(tab_pola), value = TRUE)
  tab_pola[[paste0(p, "_allele_count")]] <- rowSums(tab_pola[, cols])
}

for (p in det.pop) {
  allele_col <- paste0(p, "_allele_count")
  statut_col <- paste0(p, "_statut")
  nb_allele <- ifelse(p == "SP", 50, ifelse(p == "WP", 36, 20))
  tab_pola[[statut_col]] <- ifelse(tab_pola[[allele_col]] == 0, "Homozygous Ancestral",
                            ifelse(tab_pola[[allele_col]] == nb_allele, "Homozygous Derived", "Heterozygous"))
}

# Filter fully ancestral or derived sites
allele_count_cols <- paste0(det.pop, "_allele_count")
tab_pola$total_allele_count <- rowSums(tab_pola[, allele_count_cols])
tab_pola <- tab_pola[tab_pola$total_allele_count < 106 & tab_pola$total_allele_count > 0, ]

cat(paste0("Final number of SNPs: ", nrow(tab_pola), "\n"),
    file = "statistiques_load.txt", append = TRUE)

################################################
# Functional annotation and effect prediction
################################################

AA1_codes <- read.table("AA_codes.txt", header = TRUE)
AA2_codes <- AA1_codes
names(AA1_codes) <- c("AA1_1letter", "AA1", "AA1_Miyata_1979", "AA1_Hanada_2007_A")
names(AA2_codes) <- c("AA2_1letter", "AA2", "AA2_Miyata_1979", "AA2_Hanada_2007_A")

tab_pola <- left_join(tab_pola, AA1_codes)
tab_pola <- left_join(tab_pola, AA2_codes)

rm(AA1_codes, AA2_codes)

tab_pola$cat <- "other"
tab_pola[tab_pola$intergenic_position_status == "INTERGENIC", ]$cat <- "IG"
tab_pola[tab_pola$SNPEFF_3 == "LOW" & tab_pola$exon_position_status == "EXON", ]$cat <- "SYN"
tab_pola[tab_pola$SNPEFF_3 == "MODERATE", ]$cat <- "MIS"
tab_pola[tab_pola$SNPEFF_3 == "HIGH" & tab_pola$exon_position_status == "EXON", ]$cat <- "LOF"
tab_pola <- tab_pola[tab_pola$cat != "other", ]

cat_counts <- as.vector(table(tab_pola$cat))
cat(paste0("Number of polarized SNPs in exonic or intergenic regions: ", nrow(tab_pola), "\n"),
    file = "statistiques_load.txt", append = TRUE)

################################################
# Frequency computation and block assignment
################################################

for (p in det.pop) {
  allele_col <- paste0(p, "_allele_count")
  freq_col <- paste0(p, "_freq_der")
  nb_allele <- ifelse(p == "SP", 50, ifelse(p == "WP", 36, 20))
  tab_pola[[freq_col]] <- tab_pola[[allele_col]] / nb_allele
}

chromosome_numbers <- as.numeric(gsub("chr", "", tab_pola$Chromosome))
tab_pola <- tab_pola[order(chromosome_numbers), ]
tab_pola$block_label <- 0

JK_nbblocs <- 100
block_size <- length(tab_pola$Chromosome) / JK_nbblocs
tab_pola$block_label <- gl(JK_nbblocs, ceiling(block_size))[1:length(tab_pola$SNPbi)]

save(file = "table_pola.RData", tab_pola)
