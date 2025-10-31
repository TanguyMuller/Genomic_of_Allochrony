#!/usr/bin/env Rscript

# Load necessary libraries
library(poolfstat)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
vcf_base <- args[1]  
vcf_path <- paste0("/your/path/", vcf_base, ".vcf.gz")

# Charger les échantillons
sample.list <- as.character(read.table("/path/to/list/pool_list.list")[,1])
# Vérifier si la liste est bien chargée
if (length(sample.list) == 0) {
  stop("ERREUR : La liste des pools est vide.")
}

# Définir les tailles de pool
poolsizes <- c(52, 80, 80, 80, 80, 80, 80, 62, 78, 100, 80, 80, 68)
# Vérifier la cohérence entre poolsizes et sample.list
if (length(poolsizes) != length(sample.list)) {
  stop("ERREUR : Nombre de poolsizes différent du nombre d'échantillons.")
}

# Appliquer vcf2pooldata
tmp <- vcf2pooldata(
  vcf.file = vcf_path,
  poolsizes = poolsizes,
  poolnames = sample.list,
  min.cov.per.pool = 10
)


output_name <- gsub("_pool", "", vcf_base)
output_name <- gsub("\\.vcf\\.gz$", "", output_name)

# Créer le dataframe d'informations
info <- as.data.frame(cbind(V1 = tmp@snp.info$Chromosome, V2 = tmp@snp.info$Position))

# Écrire le fichier de sortie avec le nom correct
write.table(info[,1:2], 
            paste0("/work/project/loadexp/analyses_tanguy/Exp/list/positions_", output_name, "_10covperpool.txt"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Liste des noms de pools (correspond aux colonnes dans le VCF)
pops <- c("ES_COR", "FR_BEA", "FR_BOR", "FR_BRE", "FR_CHA", "FR_ERQ", "FR_FON", "FR_LES", "FR_MON", "FR_ORL", "FR_QUI", "FR_VEN", "IT_TRE")

# Boucle sur chaque pool
for (i in seq_along(pops)) {
  pop <- pops[i]
  # Calculer la fréquence pour le pool actuel
  freq <- as.data.frame(tmp@refallele.readcount[,i] / tmp@readcoverage[,i])
  colnames(freq) <- "freqAR"
  
  # Écrire la fréquence dans un fichier
  output_file <- paste0("freq_alleliq/frequence_", pop, "_masked_HWE_10covperpool.txt")
  write.table(freq$freqAR, output_file, row.names = FALSE, col.names = FALSE)
}

