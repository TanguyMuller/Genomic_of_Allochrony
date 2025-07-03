#!/usr/bin/env Rscript
# Script for Principal Component Analysis (PCA) on pool-seq and individual data
# Input: LD-pruned VCF files from previous filtering steps

# Load necessary libraries
library(poolfstat)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]
chr <- args[2]
variable_1 <- args[3]
variable_1 <- as.numeric(variable_1)
variable_2 <- args[4]
variable_2 <- as.numeric(variable_2)

# Load the list of samples
sample.list <- as.character(read.table("/path/to/list/all.portugal.list")[,1])

# Define the pool sizes (assuming 2 individuals per pool for diploid data)
poolsizes <- rep(2, length(sample.list))

# Apply the vcf2pooldata function to the VCF file 
tmp <- vcf2pooldata(
    vcf.file = vcf_file,
    poolsizes = poolsizes,                   
    poolnames = sample.list)

# Perform PCA on the pool and individual data
A.randombase <- randomallele.pca(tmp)

# Plot the percentage of variance explained by each principal component
pdf(paste0("Percentage_variable_explained_PC",variable_1,"~PC",variable_2,"_",chr,".pdf"), width = 600/72,height = 600/72)
perc.var <- A.randombase$perc.var
plot(perc.var, xlab = "PC", ylab = "% variance explained", type = "h", las = 2)
dev.off()

# Create a data frame with the PCA coordinates and additional information
coord.sample <- as.data.frame(A.randombase$pop.loadings)
coord.sample$pop <- as.character(read.table("/path/to/list/all.portugal.list")[,2])
coord.sample$type <- as.character(read.table("/path/to/list/all.portugal.list")[,3]) # Males, Females or Pools

# Plot the PCA results for specified PCs
library(rlang)

# Create PDF file
pdf(paste0("PCA_PC", variable_1, "~PC", variable_2, "_", chr, ".pdf"), width = 1200/72, height = 1200/72)

# Get variance percentages
pc1_var <- A.randombase$perc.var[variable_1]
pc2_var <- A.randombase$perc.var[variable_2]

# Create X and Y labels with fallback if NA
x_label <- if (!is.na(pc1_var)) {
  paste0("PC", variable_1, " (", round(pc1_var, 2), "%)")
} else {
  paste0("PC", variable_1)
}
y_label <- if (!is.na(pc2_var)) {
  paste0("PC", variable_2, " (", round(pc2_var, 2), "%)")
} else {
  paste0("PC", variable_2)
}

# Create the plot
p <- ggplot(coord.sample, aes(x = !!sym(paste0("V", variable_1)),
                                y = !!sym(paste0("V", variable_2)),
                                color = pop, shape=type)) +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dotted") +
  geom_vline(xintercept = 0, color = "black", size = 0.3, linetype = "dotted") +
  geom_point(size = 4) +
  scale_shape_manual(values = c(24, 19, 8)) + 
  scale_color_manual(values = c("brown1", "cadetblue3", "chartreuse3", "darkolivegreen4", "darkgreen", "darkorchid", 
                                "darkorange4", "darkorange2", "goldenrod2", "tomato3")) + 
  xlab(x_label) +
  ylab(y_label) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.text = element_text(color = "black", size = 22),
    axis.title = element_text(size = 22),
    axis.title.y = element_text(angle = 90),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.position = "none",
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 22)
  )

# Print the plot to PDF
print(p)

# Close PDF file
dev.off()

###Commands
# Rscript 5_PCA.R all.portugal.freebayes_filt_maf005_ldprune.vcf.gz Auto 1 2
# Rscript 5_PCA.R all.portugal.freebayes_filt_maf005_ldprune.vcf.gz chrZ 1 2
