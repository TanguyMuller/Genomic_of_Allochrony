# Script for perform two way ANOVA for test for π difference between 10 populations and between autosomes and the Z chromosome

# Set working directory
setwd("/your/way/")

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)
library(patchwork)
library(remotes)

# ============================================================
# 1. Load and combine π data from autosomes and chromosome Z
# ============================================================

# Generate the list of chromosome names
autosomes <- paste0("chr", 1:49)
chrZ <- "chrZ"

# Define the directories for autosomes and chromosome Z
autosome_dirs <- paste0("/your/way/", autosomes)
chrZ_dir <- paste0("/your/way/", chrZ)

# Initialize empty data frames to store results
combined_data_autosomes <- data.frame()
combined_data_chrZ <- data.frame()

# ---- Read and merge all autosome files ----
for (dir in autosome_dirs) {
  file_path <- file.path(dir, "pixy_pi.txt")
  pi <- read.table(file_path, header = TRUE)
  combined_data_autosomes <- rbind(combined_data_autosomes, pi)
}

# ---- Read and merge chromosome Z file ----
for (dir in chrZ_dir) {
  file_path <- file.path(dir, "pixy_pi.txt")
  pi <- read.table(file_path, header = TRUE)
  combined_data_chrZ <- rbind(combined_data_chrZ, pi)
}

# ============================================================
# 2. Combine autosome and chrZ data into one dataframe
# ============================================================

data_all <- rbind(combined_data_autosomes, combined_data_chrZ)

# Assign chromosome labels
data_all$chr <- c(rep("Auto", 58180), rep("chrZ", 2900))

# Define population order for consistent plotting
desired_order <- c("LSP", "WP", "FU", "VI", "VA", "CA", "TA", "GR", "F1", "LateLSP")
data_all$pop <- factor(data_all$pop, levels = desired_order)
data_all$pop_chr <- paste(data_all$pop, data_all$chr, sep = "_")

# ============================================================
# 3. Two-way ANOVA (population × chromosome)
# ============================================================

two_way <- data_all %>% anova_test(avg_pi ~ pop * chr)
two_way

# Post-hoc pairwise comparisons with Bonferroni correction
pwc <- data_all %>%
  emmeans_test(avg_pi ~ pop_chr, p.adjust.method = "bonferroni")
pwc

# ============================================================
# 4. Ratio of π between chrZ and autosomes
# ============================================================

set.seed(2)

# Randomly sample equal number of windows (n=290) per population
df_ratio_auto <- combined_data_autosomes %>%
  group_by(pop) %>%
  sample_n(290) %>%
  select(c(pop, avg_pi)) %>%
  ungroup()

df_ratio_Z <- combined_data_chrZ %>%
  group_by(pop) %>%
  sample_n(290) %>%
  select(c(pop, avg_pi)) %>%
  ungroup()

# Ensure both dataframes are in the same population order
df_ratio_auto <- df_ratio_auto[order(df_ratio_auto$pop), ]
df_ratio_Z <- df_ratio_Z[order(df_ratio_Z$pop), ]

# Combine and compute the ratio (chrZ / autosome)
df_ratio_combined <- bind_cols(df_ratio_auto, df_ratio_Z) %>%
  mutate(ratio_avg_pi = avg_pi...4 / avg_pi...2) %>%
  select(pop...1, ratio_avg_pi)

colnames(df_ratio_combined)[1] <- "pop"

# Compute median ratio per population
medians <- df_ratio_combined %>%
  group_by(pop) %>%
  summarize(median_ratio = median(ratio_avg_pi))

# Reorder populations by their median ratios
df_ratio_combined$pop <- factor(df_ratio_combined$pop, 
                                levels = medians$pop[order(medians$median_ratio)])

# ============================================================
# 5. One-way ANOVA on Pi ratios
# ============================================================

one_way <- df_ratio_combined %>% anova_test(ratio_avg_pi ~ pop)
one_way

# Post-hoc pairwise test with Bonferroni correction
pwc <- df_ratio_combined %>%
  emmeans_test(ratio_avg_pi ~ pop, p.adjust.method = "bonferroni")
pwc

# Save the full workspace
save.image("PI.RData")

# ============================================================
# 6. Plot 1: Average π per population (Autosomes vs chrZ)
# ============================================================

p1 <- ggplot(data_all, aes(reorder(interaction(pop, chr), avg_pi), y = avg_pi, color = pop)) +
  geom_boxplot(outliers = FALSE) +
  scale_color_manual(values = c("brown1", "cadetblue3", "chartreuse3", "darkolivegreen4", 
                                "darkgreen", "darkorchid", "darkorange4", "darkorange2", 
                                "goldenrod2", "tomato3")) +

  # ---- Manual annotation of population names ----
  # You may adjust Y positions to avoid overlap depending on your data
  annotate("text", x = 1, y = -0.0003, label = "SP", size = 7, hjust = 0.5) +
  annotate("text", x = 2, y = -0.0003, label = "LSP", size = 7, hjust = 0.5) +
  annotate("text", x = 3, y = -0.0003, label = "CA", size = 7, hjust = 0.5) +
  annotate("text", x = 4, y = -0.0003, label = "WP", size = 7, hjust = 0.5) +
  annotate("text", x = 5, y = -0.0003, label = "F1", size = 7, hjust = 0.5) +
  annotate("text", x = 6, y = -0.0003, label = "GR", size = 7, hjust = 0.5) +
  annotate("text", x = 7, y = -0.0003, label = "VI", size = 7, hjust = 0.5) +
  annotate("text", x = 8, y = -0.0003, label = "VA", size = 7, hjust = 0.5) +
  annotate("text", x = 9, y = -0.0003, label = "TA", size = 7, hjust = 0.5) +
  annotate("text", x = 10, y = -0.0003, label = "FU", size = 7, hjust = 0.5) +

  # ---- Segment labels for chrZ vs Autosomes ----
  annotate("text", x = 5, y = -0.00085, label = "chrZ", size = 7, hjust = 0.5) +
  annotate("text", x = 15, y = -0.00085, label = "Auto", size = 7, hjust = 0.5) +
  geom_segment(aes(x = 0.75, xend = 10.25, y = -0.0006, yend = -0.0006), color = "black") +
  geom_segment(aes(x = 10.75, xend = 20.25, y = -0.0006, yend = -0.0006), color = "black") +

  # ---- Final formatting ----
  labs(x = NULL, y = expression(pi), title = bquote("Average " ~ pi)) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 17),
    plot.title = element_text(hjust = 0.5, size = 17, face = "bold"),
    legend.position = "none"
  )

# ============================================================
# 7. Plot 2: Ratio of chrZ/Autosome π per population
# ============================================================

p2 <- ggplot(df_ratio_combined, aes(x = pop, y = ratio_avg_pi, color = pop)) +
  geom_boxplot(outliers = FALSE) +
  scale_color_manual(values = c("brown1", "tomato3", "darkorchid", "cadetblue3", "darkgreen",
                                "darkolivegreen4", "darkorange2", "goldenrod2", 
                                "chartreuse3", "darkorange4")) +
  
  labs(x = NULL, 
       y = bquote("Ratio of chrZ to autosomal median" ~ pi), 
       title = bquote("Average Ratio of chrZ to Autosomal Median" ~ pi)) +
  
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 17, face = "bold"),
    legend.position = "none"
  )

# ============================================================
# 8. Combine plots and save
# ============================================================

combined_plot <- p1 + p2 + plot_layout(ncol = 2, widths = c(2, 1))

# Save high-resolution figure
ggsave("combined_plot.png", plot = combined_plot, width = 20, height = 8, units = "in", dpi = 600)
