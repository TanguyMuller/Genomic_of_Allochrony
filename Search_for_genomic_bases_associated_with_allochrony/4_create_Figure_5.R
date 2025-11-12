#!/usr/bin/env Rscript

# Script to generate all panels (A–D) of Figure 5 from BayPass and population genetics results.

library(poolfstat)
library(tidyr)
library(dplyr)
library(ggpubr)
library(ape)
library(igraph)
library(ggplot2)
library(zoo)
library(svglite)
library(scales)
library(writexl)
library(openxlsx)

# ================================================================
# 1. Panel A – Genome-wide local score plot
# ================================================================

# Load BayPass combined results and helper functions
load("baypass_combine.RData")
source("baypass_utils.R")

# Load circadian genes and compute their midpoints
genes <- read.xlsx("genes_circadian.xlsx")
genes$mid <- (genes$start + genes$end) / 2

# Create output PDF file for Figure 5
pdf("Figure_5.pdf", width = 2000/72, height = 2100/72)
layout(matrix(c(1,2,3,4,5,6,7,8), nrow = 4, byrow = TRUE),
       widths = c(0.04, 0.96), heights = c(1,1,1,1))
par(oma = c(5, 5, 5, 5))

# -------------------------------
# Panel A label
# -------------------------------
par(mar = c(6, 1, 3, 0))
plot(0, 0, type = "n", axes = FALSE)
text(0.5, 0.95, "A)", cex = 5)

# Compute local scores from BayPass outputs
par(mar = c(6, 10, 3, 4))
t <- compute.local.scores(
  snp.position = df_combine[,1:2],
  snp.pi = df_pi$M_P,
  snp.pvalue = df_combine$log,
  xi = 2,
  min.maf = 0.2,
  plot.pvalhist = FALSE,
  manplot = FALSE
)

chrInf <- t$chr_info
baypass_result <- t$res.local.scores
tmp_chr <- unique(baypass_result$chr)

# Define alternating colors for chromosomes (grey/black; red for Z)
col_chr <- c(rep(c("grey", "black"), 24), "grey", "red3")

# Filter data to reduce file size (keep significant + random subset)
set.seed(42)
baypass_filtered <- baypass_result
chrZ_rows <- baypass_filtered$chr == "Z"
baypass_filtered <- baypass_filtered[
  chrZ_rows & (baypass_filtered$lindley > 0.8 | runif(sum(chrZ_rows)) < 0.005) |
  !chrZ_rows & (baypass_filtered$lindley > 0.5 | runif(sum(!chrZ_rows)) < 0.01),
]

# Compute chromosome midpoints for x-axis labels
pos_chr <- rep(0, length(tmp_chr))
tmp_nmrk <- table(baypass_filtered$chr)[tmp_chr]
pos_mrk <- cumsum(tmp_nmrk)
pos_chr[1] <- floor(pos_mrk[1]/2)
if (length(tmp_chr) > 1) {
  for (i in 2:length(tmp_chr)) {
    pos_chr[i] <- pos_mrk[i-1] + floor(tmp_nmrk[i]/2)
  }
}

# Plot Lindley local scores per SNP
data <- baypass_filtered$lindley
plot(data, pch = 16, las = 1, ylim = c(0, 85), xlim = c(0, length(data)),
     col = col_chr[baypass_filtered$chr], xaxt = "n",
     ylab = expression(C[2] ~ "Local-score"), cex.axis = 3)

# Add chromosome labels
axis(1, at = pos_chr, labels = FALSE)
selected_labels <- tmp_chr[seq(1, length(tmp_chr), by = 2)]
axis(1, at = pos_chr[seq(1, length(tmp_chr), by = 2)], 
     labels = selected_labels, las = 1, cex.axis = 3, tck = -0.02)
mtext("Chromosome", side = 1, line = 8, cex = 2.5)

# Highlight circadian genes
for (i in seq_len(nrow(genes))) {
  idx <- which(baypass_filtered$chr == genes$chr[i])
  if (length(idx) == 0) next
  mid_pos <- genes$mid[i]
  points(which.min(abs(baypass_filtered$pos[idx] - mid_pos)),
         max(baypass_filtered$lindley[idx]),
         col = alpha("black", 0.6), pch = 8, cex = 1.5)
}

# Add chromosome-specific significance thresholds
for (i in 1:nrow(chrInf)) {
  lines(x = range(which(baypass_filtered$chr == chrInf$chr[i])),
        y = rep(chrInf$th01[i], 2), col = "orange", lty = 2)
}

# ================================================================
# 2. Panel B – Local score plot for chromosome Z
# ================================================================

# Panel B label
par(mar = c(3, 1, 6, 0), mgp = c(6, 1, 0))
plot(0, 0, type = "n", axes = FALSE)
text(0.5, 0.95, "B)", cex = 5, font = 1)

# Load BayPass results for chromosome Z
par(mar = c(3, 10, 6, 4), mgp = c(6, 1, 0))

# Load circadian genes for chromosome Z
genesZ <- read.xlsx("genes_circadian_Z.xlsx")
genesZ$mid <- (genesZ$start + genesZ$end) / 2
matching_lines <- genesZ$mid  

# Measure the local-score on the Z chromosome
local_score_fb_c2_SP_WP_Z<-compute.local.scores(snp.position = Z[,1:2],   
                        snp.pi = Z_pi$M_P,                     
                        snp.pvalue=Z$log,            
                        xi=2,                        
                        min.maf=0.2,plot.pvalhist=F,        
                        manplot=F)

q1 <- quantile(local_score_fb_c2_SP_WP_Z$res.local.scores$lindley, 0.99)

# Filter to lighten the plot
set.seed(42)
local_score_fb_c2_SP_WP_Z$res.local.scores <- local_score_fb_c2_SP_WP_Z$res.local.scores[
  local_score_fb_c2_SP_WP_Z$res.local.scores$lindley > 1 | 
    runif(nrow(local_score_fb_c2_SP_WP_Z$res.local.scores)) < 0.08, 
]

# Base empty plot
plot(1,1, pch = 16, las = 1, cex = 1.2,
     ylim = c(0,85), xlim = c(0,30000000),
     col = "white", xaxt = "n",
     ylab = expression(C[2] ~ "Local-score"), cex.axis = 3, cex.lab = 3)

# Add x-axis ticks
pos_chr <- seq(0, 30000000, by = 5000000)
tmp_chr <- as.character(pos_chr)
axis(1, at = pos_chr, labels = tmp_chr, las = 1, cex.axis = 3, tck = -0.02, padj = 1)

# Gene labels (positions manually adjusted for clarity)
# Colors: black = circadian, purple = non-circadian
for (i in seq_along(matching_lines)) {
  mid <- matching_lines[i]
  color <- if (genesZ$circadian[i] == "yes") alpha("black", 0.5) else alpha("purple", 0.5)
  segments(mid, genesZ$start[i], mid, genesZ$end[i], col = color, lwd = 2)
  text(mid, genesZ$end[i] + 5, labels = genesZ$name[i], srt = 30, adj = 1, cex = 2.5, col = color)
}

# Plot local scores and significance threshold
points(local_score_fb_c2_SP_WP_Z$res.local.scores$pos,
       local_score_fb_c2_SP_WP_Z$res.local.scores$lindley, pch = 16, col = "red3")
lines(x = range(local_score_fb_c2_SP_WP_Z$res.local.scores$pos), y = rep(q1, 2),
      lty = 2, col = "orange")

# ================================================================
# 3. Panel C – FST sliding window plot (overview)
# ================================================================

# Panel C label
par(mar = c(3, 1, 3, 0))
plot(0, 0, type = "n", axes = FALSE)
text(0.1, 0.95, "C)", cex = 5, font = 1)

# Load sliding window FST data
par(mar = c(3, 10, 3, 4))
load("fst_sliding_window.RData")

# Base FST plot
plot(fst1$CumMidPos, fst1$MultiLocusFst,
     pch = 16, col = alpha("grey", 0), cex = 1,
     xaxt = "n", yaxt = "n",
     ylab = expression(F[ST] ~ "values"),
     ylim = c(0, 1.06), xlim = c(0, 30000000),
     cex.axis = 3, cex.lab = 3)

# Highlight regions of interest
rect(12000000, -0.1, 14600000, 1.2, col = rgb(0.5, 0.5, 0.5, 0.25), border = rgb(0, 0, 0, 0.35), lty = "dashed")
rect(16600000, -0.1, 18800000, 1.2, col = rgb(0.5, 0.5, 0.5, 0.25), border = rgb(0, 0, 0, 0.35), lty = "dashed")

# Add axes
axis(1, at = pos_chr, labels = tmp_chr, las = 1, cex.axis = 3, tck = -0.02)
fst_values <- seq(0, 1, 0.25)
axis(2, at = fst_values, labels = fst_values, cex.axis = 3, las = 1)

# Add smoothed FST trends
lines(rollmean(zoo(fst1$MultiLocusFst, fst_no_na1$CumMidPos), k = 200, na.rm = TRUE), col = alpha("black", 0.6), lwd = 2.5)
lines(rollmean(zoo(fst2$MultiLocusFst, fst_no_na2$CumMidPos), k = 200, na.rm = TRUE), col = alpha("red", 0.6), lwd = 2.5)
lines(rollmean(zoo(fst3$MultiLocusFst, fst_no_na3$CumMidPos), k = 200, na.rm = TRUE), col = alpha("green", 0.6), lwd = 2.5)
lines(rollmean(zoo(fst4$MultiLocusFst, fst_no_na4$CumMidPos), k = 200, na.rm = TRUE), col = alpha("orange", 0.6), lwd = 2.5)
lines(rollmean(zoo(fst5$MultiLocusFst, fst_no_na5$CumMidPos), k = 200, na.rm = TRUE), col = alpha("purple", 0.6), lwd = 2.5)

# Legend
legend("topleft",
       legend = c("LSP~LWP", "LSP~LateLSP","LWP~LateLSP", "LSP~FU", "LWP~FU"),
       col = c("black", "red", "green", "orange", "purple"),
       lwd = 1.5, bty = "n", cex = 3, ncol = 2)

# ================================================================
# 4. Panel D – FST detail plot (chromosome Z)
# ================================================================

# Panel D label
par(mar = c(7, 1, 3, 0))
plot(0, 0, type = "n", axes = FALSE)
text(0.5, 0.95, "D)", cex = 5, font = 1)

# Plot detailed FST values
par(mar = c(7, 10, 3, 4))
plot(fst1$CumMidPos, fst1$MultiLocusFst,
     pch = 16, col = rgb(0.5, 0.5, 0.5, 0.25), cex = 1.5,
     xaxt = "n", yaxt = "n",
     ylab = expression(F[ST] ~ "values"),
     xlim = c(0, 30000000),
     cex.axis = 3, cex.lab = 3)

# Add axes
axis(1, at = pos_chr, labels = tmp_chr, las = 1, cex.axis = 3, tck = -0.02)
axis(2, at = fst_values, labels = fst_values, cex.axis = 3, las = 1)

# Add smoothed FST trend
lines(rollmean(zoo(fst1$MultiLocusFst, order.by = fst1$CumMidPos),
               k = 200, na.rm = TRUE),
      col = alpha("black", 0.6), lwd = 1.5)

# Common X-axis label
mtext("Position on chromosome Z (bp)", side = 1, line = 8, cex = 2.5)

# Close the PDF device
dev.off()
