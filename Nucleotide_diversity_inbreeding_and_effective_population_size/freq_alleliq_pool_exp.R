#!/usr/bin/env Rscript

# Load necessary libraries
library(poolfstat)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
vcf_base <- args[1]  
vcf_path <- paste0("/your/path/", vcf_base, ".vcf.gz")

# Load samples
sample.list <- as.character(read.table("/path/to/list/pool.portugal.list")[,1])

# Define size pool
poolsizes <- c(100,100,80,40,52,70,74,70,34,32)

# Filt sites with fewer than 10 reads per pool
tmp <- vcf2pooldata(
  vcf.file = vcf_path,
  poolsizes = poolsizes,
  poolnames = sample.list,
  min.cov.per.pool = 10
)

# Create data frame for positions informations
info <- as.data.frame(cbind(V1 = tmp@snp.info$Chromosome, V2 = tmp@snp.info$Position))
write.table(info[,1:2], 
            paste0("/your/path/to/list/positions_", vcf_base, "_10covperpool.txt"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Estimate alleliq frequencies on each pool
for (i in seq_along(sample.list)) {
  pop <- pops[i]
  freq <- as.data.frame(tmp@refallele.readcount[,i] / tmp@readcoverage[,i])
  colnames(freq) <- "freqAR"
  output_file <- paste0("freq_alleliq/frequence_", pop, "_masked_HWE_10covperpool.txt")
  write.table(freq$freqAR, output_file, row.names = FALSE, col.names = FALSE)
}

