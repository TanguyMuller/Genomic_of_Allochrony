# Script to perform genetic load analysis (proportion of masked and realized load, Rxy, SFS and Joint SFS) 
# based on polarized SNP dataset

library(tidyverse)
library(dplyr)
library(openxlsx)
library(writexl)
library(scales)
library(patchwork)

load("table_pola.RData")
#load("Z_table_pola.RData")

# Define populations and sample sizes
det.pop <- c("SP", "WP", "FU")
nSP <- 25
#nSP <- 15
nWP <- 18
#nWP <- 12
nFU <- 10

#-----------------------------------------------------------
# 1. Compute frequency and proportion of masked and realized load
#-----------------------------------------------------------
categories <- c("IG", "SYN", "MIS", "LOF")

# Function to count genotypes per category
count_categories <- function(df, prefix, n_columns, categories) {
  results <- list()
  for (cat in categories) {
    df_cat <- subset(df, df$cat == cat)
    results_cat <- list()
    for (i in 1:n_columns) {
      col_name <- paste0(prefix, "_", i)
      freq <- table(df_cat[[col_name]])
      freq_complete <- c("0/0" = 0, "0/1" = 0, "1/1" = 0)
      freq_complete[names(freq)] <- freq
      results_cat[[col_name]] <- freq_complete
    }
    results[[cat]] <- results_cat
  }
  return(results)
}

# Count genotypes for each population
results_SP <- count_categories(tab_pola, "SP", nSP, categories)
results_WP <- count_categories(tab_pola, "WP", nWP, categories)
results_FU <- count_categories(tab_pola, "FU", nFU, categories)

results <- c(results_SP, results_WP, results_FU)

# Create dataframe with genotype counts
pop <- c(rep("LSP", nSP*12), rep("LWP", nWP*12), rep("FU", nFU*12))
state <- rep(c("Homozygous Ancestral", "Heterozygous", "Homozygous Derived"), times = ((nSP*12 + nWP*12 + nFU*12)/3)
values <- unlist(results)
individus <- c(rep(rep(1:nSP, each = 3), 4), rep(rep(1:nWP, each = 3), 4), rep(rep(1:nFU, each = 3), 4))
categorie <- c(rep("IG", (3 * nSP)), rep("SYN", (3 * nSP)), rep("MIS", (3 * nSP)), rep("LOF", (3 * nSP)),
               rep("IG", (3 * nWP)), rep("SYN", (3 * nWP)), rep("MIS", (3 * nWP)), rep("LOF", (3 * nWP)),
               rep("IG", (3 * nFU)), rep("SYN", (3 * nFU)), rep("MIS", (3 * nFU)), rep("LOF", (3 * nFU)))

nb_cat <- data.frame(pop = pop, individus = individus, state = state, values = values, categorie = categorie)

# Calculate median and proportions
df_median <- nb_cat %>%
  group_by(pop, state, categorie) %>%
  summarize(median_value = median(values, na.rm = TRUE), .groups = "drop")

total_by_pop_cat <- df_median %>%
  group_by(pop, categorie) %>%
  summarize(total = sum(median_value, na.rm = TRUE), .groups = "drop")

df_median <- df_median %>%
  left_join(total_by_pop_cat, by = c("pop", "categorie")) %>%
  mutate(proportion = median_value / total)

df_median$pop <- factor(df_median$pop, levels = c("LSP", "LWP", "FU"))
df_median$state <- factor(df_median$state, levels = c("Homozygous Ancestral", "Heterozygous", "Homozygous Derived"))

# Add total SNP counts
total_snps <- c(IG = 2348849, SYN = 47248, MIS = 29437, LOF = 374)
#total_snps <- c(IG=52430, SYN=1742, MIS=1155, LOF=12)
df_median$Category_label <- total_snps[as.character(df_median$categorie)]
df_median$categorie <- factor(df_median$categorie, levels = c("IG", "SYN", "MIS", "LOF"))

# Create line dataframe for plot
population_names <- unique(df_median$pop)
first_population <- population_names[2]

line_df_median_count <- data.frame(
  x = rep(0, length(population_names)),
  xend = rep(0, length(population_names)),
  y = ifelse(population_names == first_population, -0.002, NA),
  yend = ifelse(population_names == first_population, 1.002, NA),
  pop = population_names
)

# Plot genotype proportions
count <- ggplot(df_median, aes(x = categorie, y = proportion, fill = state)) +
  geom_bar(stat = "identity", color = "black") +
  geom_segment(data = line_df_median_count, 
               aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE,
               color = "black", size = 1) +
  geom_text(data = df_median,
            aes(x = categorie, y = 1.02, label = Category_label),
            inherit.aes = FALSE,
            size = 3,
            vjust = 0, color = "grey20") +
  facet_wrap(~ pop) +
  scale_fill_manual(
    values = c(
      "Homozygous Ancestral" = "#999999",
      "Heterozygous" = "#E69F00",
      "Homozygous Derived" = "#56B4E9"
    )
  ) +
  scale_y_continuous(limits = c(-0.01, 1.05), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 14, color = "grey20"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.spacing.x = unit(1.5, "cm"),
    legend.text = element_text(color = "grey20", size = 10),
    legend.title = element_blank()
  ) +
  labs(x = NULL, y = "Proportion", fill = NULL)

count
save(count, df_median, line_df_median_count, file = "count.RData")
#save(count, df_median, line_df_median_count, file = "Z_count.RData")

# Total derived allele count
df_new <- nb_cat %>%
  filter(state %in% c("Heterozygous", "Homozygous Derived")) %>%
  mutate(weighted_count = case_when(
    state == "Heterozygous" ~ values * 1,
    state == "Homozygous Derived" ~ values * 2
  )) %>%
  group_by(pop, individus, categorie) %>%
  summarise(total_derived_alleles = sum(weighted_count), .groups = "drop")

#-----------------------------------------------------------
# 2. Rxy (Relative load)
#-----------------------------------------------------------
Rxy <- data.frame("IG" = 0, "SYN" = 0, "MIS" = 0, "LOF" = 0, "Data" = 0)
categories <- c("IG", "SYN", "MIS", "LOF")
JK_nbblocs <- 100

# Rxy SP vs WP
d <- 1
for (g in 1:JK_nbblocs) {
  tab_pola_g <- tab_pola[tab_pola$block_label != g & !is.na(tab_pola$cat), ]
  
  for (cat in categories) {
    tab_cat <- tab_pola_g[tab_pola_g$cat == cat, ]
    col1 <- "SP_freq_der"
    col2 <- "WP_freq_der"
    data <- tab_cat[, c(col1, col2)]
    
    data$RX <- data$SP_freq_der * (1 - data$WP_freq_der)
    data$RY <- data$WP_freq_der * (1 - data$SP_freq_der)
    
    Rxy[d, cat] <- sum(data$RX, na.rm = TRUE) / sum(data$RY, na.rm = TRUE)
    Rxy[d, "Data"] <- g
  }
  d <- d + 1
}
Rxy_SP_WP <- Rxy

# Rxy SP vs FU
d <- 1
for (g in 1:JK_nbblocs) {
  tab_pola_g <- tab_pola[tab_pola$block_label != g & !is.na(tab_pola$cat), ]
  
  for (cat in categories) {
    tab_cat <- tab_pola_g[tab_pola_g$cat == cat, ]
    col1 <- "SP_freq_der"
    col2 <- "FU_freq_der"
    data <- tab_cat[, c(col1, col2)]
    
    data$RX <- data$SP_freq_der * (1 - data$FU_freq_der)
    data$RY <- data$FU_freq_der * (1 - data$SP_freq_der)
    
    Rxy[d, cat] <- sum(data$RX, na.rm = TRUE) / sum(data$RY, na.rm = TRUE)
    Rxy[d, "Data"] <- g
  }
  d <- d + 1
}
Rxy_SP_FU <- Rxy

# Rxy WP vs FU
d <- 1
for (g in 1:JK_nbblocs) {
  tab_pola_g <- tab_pola[tab_pola$block_label != g & !is.na(tab_pola$cat), ]
  
  for (cat in categories) {
    tab_cat <- tab_pola_g[tab_pola_g$cat == cat, ]
    col1 <- "WP_freq_der"
    col2 <- "FU_freq_der"
    data <- tab_cat[, c(col1, col2)]
    
    data$RX <- data$WP_freq_der * (1 - data$FU_freq_der)
    data$RY <- data$FU_freq_der * (1 - data$WP_freq_der)
    
    Rxy[d, cat] <- sum(data$RX, na.rm = TRUE) / sum(data$RY, na.rm = TRUE)
    Rxy[d, "Data"] <- g
  }
  d <- d + 1
}
Rxy_WP_FU <- Rxy

# Combine and standardize by IG
Rxy <- rbind(Rxy_SP_WP, Rxy_SP_FU, Rxy_WP_FU)
Rxy$comparaison <- c(rep("LSP/LWP", 100), rep("LSP/FU", 100), rep("LWP/FU", 100))

for (i in 1:nrow(Rxy)) {
  Rxy[i, 1:4] <- sweep(Rxy[i, 1:4], MARGIN = 1, STATS = Rxy[i, "IG"], FUN = '/')
}

save(Rxy, file = "Rxy.RData")
#save(Rxy, file = "Z_Rxy.Rdata")

# Statistical test example
Rxy_SP_FU <- Rxy[Rxy$comparaison == "LWP/FU", ]
B <- nrow(Rxy_SP_FU)
mean_val <- mean(Rxy_SP_FU$MIS)
se_jack <- sqrt((B - 1) / B * sum((Rxy_SP_FU$MIS - mean_val)^2))
z <- (mean_val - 1) / se_jack
p <- 2 * pnorm(-abs(z))

#-----------------------------------------------------------
# 3. Site Frequency Spectrum (SFS)
#-----------------------------------------------------------

# SFS for FU
nb_categories <- (0:(2 * nFU)) / (2 * nFU)
cat_list <- c("IG", "SYN", "MIS", "LOF")
results_list_FU <- list()
k <- 1

for (j in seq_along(cat_list)) {
  cat_name <- cat_list[j]
  tab_cat <- tab_pola[tab_pola$cat == cat_name & !is.na(tab_pola$cat), ]
  total_cat <- nrow(tab_cat)
  
  for (i in nb_categories) {
    for (pop in det.pop) {
      col <- paste0(pop, "_allele_count")
      somme <- sum((tab_cat[[col]] / (2 * nFU)) == i, na.rm = TRUE)
      proportion <- somme / total_cat
      
      results_list_FU[[k]] <- data.frame(
        pop = pop,
        Frequency = proportion,
        cat = cat_name,
        bin = i
      )
      k <- k + 1
    }
  }
}

tab_pola_cat_FU <- do.call(rbind, results_list_FU)

sp <- tab_pola_cat_FU[tab_pola_cat_FU$pop == "FU", ]
cat_order <- c("IG", "SYN", "MIS", "LOF")

p4 <- ggplot(sp, aes(x = bin, y = Frequency, group = cat)) +
  theme_classic() +
  geom_line(aes(color = factor(cat, levels = cat_order)), size = 2) +
  scale_color_manual(values = c("palevioletred1", "darkorange", "red", "darkred")) +
  labs(title = "FU Autosomes", x = "Derived allele frequency", y = "Frequency", color = "") +
  #labs(title = "FU chromosome Z", x = "Derived allele frequency", y = "Frequency", color = "") +
  ylim(0, 0.5) +
  #ylim(0, 0.7) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.border = element_rect(color = "grey", fill = NA, size = 1),
    axis.line = element_line(color = "grey"),
    panel.grid.major = element_line(color = alpha("grey", 0.3)),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = c(0.5, 0.95),
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box.just = "center",
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 18)
  )

# SFS for SP
nb_categories <- (0:(2 * nSP)) / (2 * nSP)
results_list_SP <- list()
k <- 1

for (j in seq_along(cat_list)) {
  cat_name <- cat_list[j]
  tab_cat <- tab_pola[tab_pola$cat == cat_name & !is.na(tab_pola$cat), ]
  total_cat <- nrow(tab_cat)
  
  for (i in nb_categories) {
    for (pop in det.pop) {
      col <- paste0(pop, "_allele_count")
      somme <- sum(tab_cat[[col]] / (2 * nSP) == i, na.rm = TRUE)
      proportion <- somme / total_cat
      
      results_list_SP[[k]] <- data.frame(
        pop = pop,
        Frequency = proportion,
        cat = cat_name,
        bin = i
      )
      k <- k + 1
    }
  }
}

tab_pola_cat <- do.call(rbind, results_list_SP)
sp <- tab_pola_cat[tab_pola_cat$pop == "SP", ]

p2 <- ggplot(sp, aes(x = bin, y = Frequency, group = cat)) +
  theme_classic() +
  geom_line(aes(color = factor(cat, levels = cat_order)), size = 2) +
  scale_color_manual(values = c("palevioletred1", "darkorange", "red", "darkred")) +
  labs(title = "LSP Autosomes", x = "Allele count", y = "Frequency", color = "") +
  ylim(0, 0.5) +
  #labs(title = "LSP Chromosome Z", x = "Allele count", y = "Frequency", color = "") +
  #ylim(0, 0.7) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.border = element_rect(color = "grey", fill = NA, size = 1),
    axis.line = element_line(color = "grey"),
    panel.grid.major = element_line(color = alpha("grey", 0.3)),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = c(0.5, 0.95),
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box.just = "center",
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 18)
  )

# SFS for WP
nb_categories <- (0:(2 * nWP)) / (2 * nWP)
results_list_WP <- list()
k <- 1

for (j in seq_along(cat_list)) {
  cat_name <- cat_list[j]
  tab_cat <- tab_pola[tab_pola$cat == cat_name & !is.na(tab_pola$cat), ]
  total_cat <- nrow(tab_cat)
  
  for (i in nb_categories) {
    for (pop in det.pop) {
      col <- paste0(pop, "_allele_count")
      somme <- sum(tab_cat[[col]] / (2 * nWP) == i, na.rm = TRUE)
      proportion <- somme / total_cat
      
      results_list_WP[[k]] <- data.frame(
        pop = pop,
        Frequency = proportion,
        cat = cat_name,
        bin = i
      )
      k <- k + 1
    }
  }
}

tab_pola_cat_WP <- do.call(rbind, results_list_WP)
sp <- tab_pola_cat_WP[tab_pola_cat_WP$pop == "WP", ]

p3 <- ggplot(sp, aes(x = bin, y = Frequency, group = cat)) +
  theme_classic() +
  geom_line(aes(color = factor(cat, levels = cat_order)), size = 2) +
  scale_color_manual(values = c("palevioletred1", "darkorange", "red", "darkred")) +
  labs(title = "LWP Autosomes", x = "Allele count", y = "Frequency", color = "") +
  ylim(0, 0.5) +
  #labs(title = "LWP Chromosome Z", x = "Allele count", y = "Frequency", color = "") +
  #ylim(0, 0.7) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.border = element_rect(color = "grey", fill = NA, size = 1),
    axis.line = element_line(color = "grey"),
    panel.grid.major = element_line(color = alpha("grey", 0.3)),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = c(0.5, 0.95),
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box.just = "center",
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 18)
  )

# Save combined SFS plot
pdf("INTERGENIC-EXON_SFS_fixed_sites.pdf", width = 20, height = 8)
#pdf("Z_INTERGENIC-EXON_SFS_fixed_sites.pdf", width = 20, height = 8)
tryCatch({
  combined_plot <- (p2 + p3 + p4) + plot_layout(ncol = 3)
  print(combined_plot)
}, error = function(e) {
  message("Error during plotting: ", e$message)
})
dev.off()

#-----------------------------------------------------------
# 4. Joint SFS (2D frequency plots)
#-----------------------------------------------------------

SP_WP <- tab_pola[, c(182:184)]
SP_FU <- tab_pola[, c(182:183,185)]
WP_FU <- tab_pola[, c(182,184:185)]
#SP_WP <- tab_pola[, c(134:136)]
#SP_FU <- tab_pola[, c(134:135,137)]
#WP_FU <- tab_pola[, c(134,136,137)]

categories <- c("IG", "SYN", "MIS", "LOF")

# Joint SFS SP vs FU
plot_list_SP_FU <- list()

for (cat in categories) {
  subset_data <- SP_FU[SP_FU$cat == cat, ]
  
  if (nrow(subset_data) > 0) {
    temp_plot <- ggplot(subset_data, aes(x = FU_freq_der, y = SP_freq_der)) +
      stat_bin_2d(bins = 10)
    
    plot_data <- ggplot_build(temp_plot)$data[[1]]
    max_count <- max(plot_data$count, na.rm = TRUE)
    max_power <- ceiling(log10(max_count))
    max_limit <- 10^max_power
    start_power <- ifelse(max_power <= 3, 0, 1)
    powers <- seq(start_power, max_power, by = 1)
    breaks <- 10^powers
    
    p <- ggplot(subset_data, aes(x = FU_freq_der, y = SP_freq_der)) +
      stat_bin_2d(bins = 10, color = "white", size = 0.25) +
      scale_fill_gradientn(
        colours = c("blue", "yellow", "red"),
        trans = "log10",
        name = "count",
        breaks = breaks,
        labels = function(x) parse(text = paste0("10^", log10(x))),
        limits = c(ifelse(start_power == 0, 1, 10), max_limit),
        oob = scales::squish
      ) +
      coord_fixed() +
      labs(x = "Allele frequency in FU", y = "Allele frequency in LSP", title = cat) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 22, hjust = 0.5),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22, color = "black", margin = margin(r = 20)),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.5),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.key.height = unit(1.25, "cm"),
        legend.key.width = unit(0.75, "cm")
      )
    
    plot_list_SP_FU[[cat]] <- p
  }
}

# Joint SFS SP vs WP
plot_list_SP_WP <- list()

for (cat in categories) {
  subset_data <- SP_WP[SP_WP$cat == cat, ]
  
  if (nrow(subset_data) > 0) {
    temp_plot <- ggplot(subset_data, aes(x = WP_freq_der, y = SP_freq_der)) +
      stat_bin_2d(bins = 10)
    
    plot_data <- ggplot_build(temp_plot)$data[[1]]
    max_count <- max(plot_data$count, na.rm = TRUE)
    max_power <- ceiling(log10(max_count))
    max_limit <- 10^max_power
    start_power <- ifelse(max_power <= 3, 0, 1)
    powers <- seq(start_power, max_power, by = 1)
    breaks <- 10^powers
    
    p <- ggplot(subset_data, aes(x = WP_freq_der, y = SP_freq_der)) +
      stat_bin_2d(bins = 10, color = "white", size = 0.25) +
      scale_fill_gradientn(
        colours = c("blue", "yellow", "red"),
        trans = "log10",
        name = "count",
        breaks = breaks,
        labels = function(x) parse(text = paste0("10^", log10(x))),
        limits = c(ifelse(start_power == 0, 1, 10), max_limit),
        oob = scales::squish
      ) +
      coord_fixed() +
      labs(x = "Allele frequency in LWP", y = "Allele frequency in LSP", title = cat) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 22, hjust = 0.5),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22, color = "black", margin = margin(r = 20)),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.5),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.key.height = unit(1.25, "cm"),
        legend.key.width = unit(0.75, "cm")
      )
    
    plot_list_SP_WP[[cat]] <- p
  }
}

# Joint SFS WP vs FU
plot_list_WP_FU <- list()

for (cat in categories) {
  subset_data <- WP_FU[WP_FU$cat == cat, ]
  
  if (nrow(subset_data) > 0) {
    temp_plot <- ggplot(subset_data, aes(x = FU_freq_der, y = WP_freq_der)) +
      stat_bin_2d(bins = 10)
    
    plot_data <- ggplot_build(temp_plot)$data[[1]]
    max_count <- max(plot_data$count, na.rm = TRUE)
    max_power <- ceiling(log10(max_count))
    max_limit <- 10^max_power
    start_power <- ifelse(max_power <= 3, 0, 1)
    powers <- seq(start_power, max_power, by = 1)
    breaks <- 10^powers
    
    p <- ggplot(subset_data, aes(x = FU_freq_der, y = WP_freq_der)) +
      stat_bin_2d(bins = 10, color = "white", size = 0.25) +
      scale_fill_gradientn(
        colours = c("blue", "yellow", "red"),
        trans = "log10",
        name = "count",
        breaks = breaks,
        labels = function(x) parse(text = paste0("10^", log10(x))),
        limits = c(ifelse(start_power == 0, 1, 10), max_limit),
        oob = scales::squish
      ) +
      coord_fixed() +
      labs(x = "Allele frequency in FU", y = "Allele frequency in LWP", title = cat) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 22, hjust = 0.5),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22, color = "black", margin = margin(r = 20)),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.5),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.key.height = unit(1.25, "cm"),
        legend.key.width = unit(0.75, "cm")
      )
    
    plot_list_WP_FU[[cat]] <- p
  }
}

# Combine all joint SFS plots
all_plots <- c(plot_list_SP_WP, plot_list_SP_FU, plot_list_WP_FU)
combined_plot <- wrap_plots(all_plots, ncol = 4)
ggsave("combined_SFS_joint.pdf", combined_plot, width = 2000/72, height = 1200/72)
#ggsave("Z_combined_SFS_joint.pdf", combined_plot, width = 2000/72, height = 1200/72)
