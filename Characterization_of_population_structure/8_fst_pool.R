#!/usr/bin/env Rscript
# Script to calculates pairwise FST between populations for autosomes and Z chromosome
# and generates phylogenetic trees and exports FST matrices

library(poolfstat)
library(data.table)
library(dplyr)
library(openxlsx)

# Get command line arguments for VCF files
args <- commandArgs(trailingOnly = TRUE)
vcf_Auto <- args[1]  # Autosome VCF file
vcf_Z <- args[2]     # Z chromosome VCF file

# Load pool names
pool.list = as.character(read.table("/path/to/list/pool.portugal.list")[,1])

# Diploid pool sizes for pools that contains males and females
diploide_size_VI <- 20
diploide_size_VA <- 26
diploide_size_F1 <- 17

# Male proportions calculated from median coverage ratios (chrZ/chr1)
# Example: 0.32/0.5 = coverage ratio indicating sex ratio in pools
m1 <- 0.32 / 0.5  # Male proportion in pool VI
m2 <- 0.21 / 0.5  # Male proportion in pool VA
m3 <- 0.26 / 0.5  # Male proportion in pool F1

# Calculate haploid sizes for chrZ (accounting for sex chromosome dosage)
# Males have 2 copies, females have 1 copy of Z chromosome
nbrVI <- m1 * diploide_size_VI * 2 + (1 - m1) * diploide_size_VI * 1
nbrVI <- round(nbrVI)
nbrVA <- m2 * diploide_size_VA * 2 + (1 - m2) * diploide_size_VA * 1
nbrVA <- round(nbrVA)
nbrF1 <- m3 * diploide_size_F1 * 2 + (1 - m3) * diploide_size_F1 * 1
nbrF1 <- round(nbrF1)

# Pool sizes for autosomes (diploid)
haploid_size_auto=c("100","100","80","40","52","70","74","70","34","32")
# Pool sizes for Z chromosome (adjusted for sex ratio)
haploid_size_Z=c("100","100","80",nbrVI,nbrVA,"70","74","70",nbrF1,"32")

# Load autosome data
tmp <- vcf2pooldata(vcf_Auto,
                    poolsizes = haploid_size_auto,
                    poolnames = pool.list)
tmp <- pooldata.subset(tmp,pool.index = c(1,10,9,2,6,7,8,5,4,3)) # Reorder pools

# Load Z chromosome data
tmp.Z <- vcf2pooldata(vcf_Z,
                      poolsizes = haploid_size_Z,
                      poolnames = pool.list)
tmp.Z <- pooldata.subset(tmp.Z,pool.index = c(1,10,9,2,6,7,8,5,4,3)) # Reorder pools

# Calculate pairwise FST for autosomes
fst.Auto <- compute.pairwiseFST(tmp,method = "Identity",nsnp.per.bjack.block = 1000)

# Calculate pairwise FST for Z chromosome
fst.Z <- compute.pairwiseFST(tmp.Z, method = "Identity", nsnp.per.bjack.block = 1000)

# Create phylogenetic trees based on FST distances
library(ape)
library(igraph)
library(ggplot2)

# Autosome tree
fst.Auto.tree <- fst.Auto@PairwiseFSTmatrix
tree <- nj(dist(fst.Auto.tree))  # Neighbor-joining tree
plot(tree,cex=1,lwd=10,"phylogram")
tree_ladderize=ladderize(tree)
plot(tree_ladderize,cex=1,lwd=10,'phylogram')
tree_root <- root(tree_ladderize, out=7)  # Root tree with outgroup
tree_root <- ladderize(tree_root)

# Save autosome tree
png("tree_dist_fst_auto.png",height = 800,width = 600)
plot(tree_root,cex=1,'phylogram')
dev.off()

# Z chromosome tree
fst.Z.tree <- fst.Z@PairwiseFSTmatrix
tree <- nj(dist(fst.Z.tree))  # Neighbor-joining tree
plot(tree,cex=1,lwd=10,"phylogram")
tree_ladderize=ladderize(tree)
plot(tree_ladderize,cex=1,lwd=10,'phylogram')
tree_root <- root(tree_ladderize, out=8)  # Root tree with outgroup
tree_root <- ladderize(tree_root)

# Save Z chromosome tree
png("tree_dist_fst_Z.png",height = 800,width = 600)
plot(tree_root,cex=1,'phylogram')
dev.off()

# Export FST matrices to Excel files
write.xlsx(fst.Auto.tree, file = "matrice_fst_auto.xlsx")
write.xlsx(fst.Z.tree, file = "matrice_fst_Z.xlsx")

###Command to run 8_fst_pool.R
#Rscript 8_fst_pool.R pool pool.portugal.freebayes_filt_maf005_ldprune.vcf.gz pool.portugal.freebayes.chrZ_filt_maf005_ldprune.vcf.gz
