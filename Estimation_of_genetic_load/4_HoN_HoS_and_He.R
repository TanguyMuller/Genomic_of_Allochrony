# Heterozygosity and HoN/HoS analysis

library(tidyverse)
library(dplyr)

nSP=25
#nSP=15
nWP=18
#nWP=12
nFU=10

# Load genotype data
all <- read.delim("vcf/load.all.GT.txt", skip=1, header=FALSE)
#all <- read.delim("vcf/Z.load.all.GT.txt", skip=1, header=FALSE)
SP <- paste0("SP_", 1:nSP)
WP <- paste0("WP_", 1:nWP)
FU <- paste0("FU_", 1:nFU)

colnames(all) <- c("Chromosome", 
                   "position", 
                   "REF", 
                   "ALT", 
                   SP, 
                   WP,
                   FU,
                   "SNPEFF_2", 
                   "SNPEFF_3", 
                   "SNPEFF_11")

all$SNPEFF_2 <- as.factor(all$SNPEFF_2)
all$SNPEFF_3 <- as.factor(all$SNPEFF_3)
all$SNPEFF_11 <- as.factor(all$SNPEFF_11)

# Add annotation columns
all$AA1 <- as.factor(str_sub(all$SNPEFF_11, 3, 5))
all$AA2 <- as.factor(ifelse(all$SNPEFF_3=="HIGH", NA, str_sub(all$SNPEFF_11, -3, -1)))

# Add genomic region status
all$exon_position_status <- as.factor(ifelse(is.na(all$SNPEFF_11)==TRUE, "not_exon", "EXON"))
all$intron_position_status <- as.factor(ifelse(all$SNPEFF_2=="intron_variant", "INTRON", "not_intron"))
all$intergenic_position_status <- as.factor(ifelse(all$SNPEFF_2=="intergenic_region", "INTERGENIC", "not_intergenic"))

all <- all[!is.na(all$intergenic_position_status), ]
all$SNPbi <- paste(all$Chromosome,"__",all$position, sep="")

# Add categories annotation
all$cat1 <- "other" 
all[!is.na(all$intergenic_position_status) & all$intergenic_position_status=="INTERGENIC" , ]$cat1 <- "IG"
all[all$SNPEFF_3=="LOW" & all$exon_position_status=="EXON",]$cat1 <- "SYN"
all[all$SNPEFF_3=="MODERATE",]$cat1 <- "MIS"
all[all$SNPEFF_3=="HIGH" & all$exon_position_status=="EXON",]$cat1 <- "LOF"
table(all$cat1)

#-----------------------------------------------------------
# 1. PiN/PiS ratio per individual (based on observed heterozygosity)
#-----------------------------------------------------------

# Calculate for SP population
results_SP <- data.frame()

for (i in 1:nSP) {
  colname <- paste0("SP_", i)
  
  # Calculate heterozygosity for each category
  NS <- all[!is.na(all$cat1) & all$cat1 == "MIS", ]
  Ho_NS <- sum(NS[[colname]] == "0/1") / nrow(NS)
  
  S <- all[!is.na(all$cat1) & all$cat1 == "SYN", ]
  Ho_S <- sum(S[[colname]] == "0/1") / nrow(S)
  
  NSH <- all[!is.na(all$cat1) & all$cat1 == "LOF", ]
  Ho_NSH <- sum(NSH[[colname]] == "0/1") / nrow(NSH)
  
  individual_results <- data.frame(
    ind = colname,
    sum_PiS_Ho = Ho_S,
    sum_PiN_Ho = Ho_NS,
    sum_PiNSH_Ho = Ho_NSH,
    PiN_PiS_Ho = ifelse(Ho_S != 0, Ho_NS / Ho_S, NA),
    PiNH_PiS_Ho = ifelse(Ho_S != 0, Ho_NSH / Ho_S, NA)
  )
  
  results_SP <- rbind(results_SP, individual_results)
}

print(results_SP)

# Calculate for WP population
results_WP <- data.frame()

for (i in 1:nWP) {
  colname <- paste0("WP_", i)
  
  NS <- all[!is.na(all$cat1) & all$cat1 == "MIS", ]
  Ho_NS <- sum(NS[[colname]] == "0/1") / nrow(NS)
  
  S <- all[!is.na(all$cat1) & all$cat1 == "SYN", ]
  Ho_S <- sum(S[[colname]] == "0/1") / nrow(S)
  
  NSH <- all[!is.na(all$cat1) & all$cat1 == "LOF", ]
  Ho_NSH <- sum(NSH[[colname]] == "0/1") / nrow(NSH)
  
  individual_results <- data.frame(
    ind = colname,
    sum_PiS_Ho = Ho_S,
    sum_PiN_Ho = Ho_NS,
    sum_PiNSH_Ho = Ho_NSH,
    PiN_PiS_Ho = ifelse(Ho_S != 0, Ho_NS / Ho_S, NA),
    PiNH_PiS_Ho = ifelse(Ho_S != 0, Ho_NSH / Ho_S, NA)
  )
  
  results_WP <- rbind(results_WP, individual_results)
}

print(results_WP)

# Calculate for FU population
results_FU <- data.frame()

for (i in 1:nFU) {
  colname <- paste0("FU_", i)
  
  NS <- all[!is.na(all$cat1) & all$cat1 == "MIS", ]
  Ho_NS <- sum(NS[[colname]] == "0/1") / nrow(NS)
  
  S <- all[!is.na(all$cat1) & all$cat1 == "SYN", ]
  Ho_S <- sum(S[[colname]] == "0/1") / nrow(S)
  
  NSH <- all[!is.na(all$cat1) & all$cat1 == "LOF", ]
  Ho_NSH <- sum(NSH[[colname]] == "0/1") / nrow(NSH)
  
  individual_results <- data.frame(
    ind = colname,
    sum_PiS_Ho = Ho_S,
    sum_PiN_Ho = Ho_NS,
    sum_PiNSH_Ho = Ho_NSH,
    PiN_PiS_Ho = ifelse(Ho_S != 0, Ho_NS / Ho_S, NA),
    PiNH_PiS_Ho = ifelse(Ho_S != 0, Ho_NSH / Ho_S, NA)
  )
  
  results_FU <- rbind(results_FU, individual_results)
}

print(results_FU)

# Combine results
piN_piS <- rbind(results_SP, results_WP, results_FU)
piN_piS$pop <- c(rep("LSP", nSP), rep("LWP", nWP), rep("FU", nFU))

piN_piS$pop <- factor(piN_piS$pop, levels = c("LSP", "LWP", "FU"))
color_pinpis <- c("brown1", "cadetblue3", "chartreuse3")

# Plot PiN/PiS ratio
pin_pis <- ggplot(data=piN_piS, aes(pop, PiN_PiS_Ho, fill=pop)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.8), 
              size = 2, color = "black", alpha = 0.5) +
  scale_fill_manual(values = color_pinpis) +
  theme(
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    axis.text.x = element_text(size = 26),
    axis.text.y = element_text(size = 26),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm")
  ) +
  labs(x = NULL, y = "NS/S Ratio")

pin_pis

save(pin_pis, piN_piS, file = "pin_pis_all.RData")
#save(pin_pis, piN_piS, file = "Z_pin_pis_all.RData")

# Statistical tests for PiN/PiS
load("pin_pis_all.RData")
#load("Z_pin_pis_all.RData")
shapiro_test_results <- by(piN_piS$PiN_PiS_Ho, interaction(piN_piS$pop), shapiro.test)
shapiro_test_results

ggplot(piN_piS, aes(sample = PiN_PiS_Ho)) +
  facet_wrap(~ pop) +
  stat_qq() +
  stat_qq_line() +
  theme_minimal() +
  labs(title = "QQ Plots by Population",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")

kruskal.test(PiN_PiS_Ho ~ pop, data = piN_piS)
library(FSA)
dunnTest(PiN_PiS_Ho ~ pop, data = piN_piS, method = "bonferroni")

#-----------------------------------------------------------
# 2. Expected heterozygosity (He) and observed heterozygosity (Ho)
#-----------------------------------------------------------

# Calculate for SP
geno_cols <- grep("^SP_", colnames(all), value = TRUE)
geno_mat <- as.matrix(all[, geno_cols])

alt_counts <- matrix(NA_integer_, nrow = nrow(geno_mat), ncol = ncol(geno_mat),
                     dimnames = dimnames(geno_mat))
gc()

# Convert genotypes to allele counts
alt_counts[geno_mat == "0/0"] <- 0
alt_counts[geno_mat == "0/1"] <- 1
alt_counts[geno_mat == "1/0"] <- 1
alt_counts[geno_mat == "1/1"] <- 2

# Calculate allele frequencies
sum_alt <- rowSums(alt_counts, na.rm = TRUE)
n_alleles <- rowSums(!is.na(alt_counts)) * 2
p <- (n_alleles - sum_alt) / n_alleles
q <- 1 - p

# Expected heterozygosity (He = 2pq)
all$pi_SP <- 2 * p * q

# Observed heterozygosity (Ho)
all$Ho_SP <- rowSums(all[, geno_cols] == "0/1", na.rm = TRUE) / length(geno_cols)

# Calculate by category
piS <- mean(all$Ho_SP[all$cat1 == "SYN"], na.rm = TRUE)
PiS <- mean(all$pi_SP[all$cat1 == "SYN"], na.rm = TRUE)
piN <- mean(all$Ho_SP[all$cat1 == "MIS"], na.rm = TRUE)
PiN <- mean(all$pi_SP[all$cat1 == "MIS"], na.rm = TRUE)
piIG <- mean(all$Ho_SP[all$cat1 == "IG"], na.rm = TRUE)
PiIG <- mean(all$pi_SP[all$cat1 == "IG"], na.rm = TRUE)
piH <- mean(all$Ho_SP[all$cat1 == "LOF"], na.rm = TRUE)
PiH <- mean(all$pi_SP[all$cat1 == "LOF"], na.rm = TRUE)
PiA <- mean(all$pi_SP, na.rm = TRUE)
PiN_PiS <- PiN / PiS
PiH_Pis <- PiH / PiS

results <- data.frame(
  pop = "SP",
  he_s = PiS,
  he_ns = PiN,
  he_h = PiH,
  he_ig = PiIG,
  he_all = PiA,
  ho_s = piS,
  ho_ns = piN,
  ho_h = piH,
  ho_ig = piIG
)

# Calculate for WP
geno_cols <- grep("^WP_", colnames(all), value = TRUE)
geno_mat <- as.matrix(all[, geno_cols])

alt_counts <- matrix(NA_integer_, nrow = nrow(geno_mat), ncol = ncol(geno_mat),
                     dimnames = dimnames(geno_mat))
gc()

alt_counts[geno_mat == "0/0"] <- 0
alt_counts[geno_mat == "0/1"] <- 1
alt_counts[geno_mat == "1/0"] <- 1
alt_counts[geno_mat == "1/1"] <- 2

sum_alt <- rowSums(alt_counts, na.rm = TRUE)
n_alleles <- rowSums(!is.na(alt_counts)) * 2
p <- (n_alleles - sum_alt) / n_alleles
q <- 1 - p

all$pi_WP <- 2 * p * q
all$Ho_WP <- rowSums(all[, geno_cols] == "0/1", na.rm = TRUE) / length(geno_cols)

piS <- mean(all$Ho_WP[all$cat1 == "SYN"], na.rm = TRUE)
PiS <- mean(all$pi_WP[all$cat1 == "SYN"], na.rm = TRUE)
piN <- mean(all$Ho_WP[all$cat1 == "MIS"], na.rm = TRUE)
PiN <- mean(all$pi_WP[all$cat1 == "MIS"], na.rm = TRUE)
piIG <- mean(all$Ho_WP[all$cat1 == "IG"], na.rm = TRUE)
PiIG <- mean(all$pi_WP[all$cat1 == "IG"], na.rm = TRUE)
piH <- mean(all$Ho_WP[all$cat1 == "LOF"], na.rm = TRUE)
PiH <- mean(all$pi_WP[all$cat1 == "LOF"], na.rm = TRUE)
PiA <- mean(all$pi_WP, na.rm = TRUE)
PiN_PiS <- PiN / PiS
PiH_Pis <- PiH / PiS

results_WP <- data.frame(
  pop = "WP",
  he_s = PiS,
  he_ns = PiN,
  he_h = PiH,
  he_ig = PiIG,
  he_all = PiA,
  ho_s = piS,
  ho_ns = piN,
  ho_h = piH,
  ho_ig = piIG
)

# Calculate for FU
geno_cols <- grep("^FU_", colnames(all), value = TRUE)
geno_mat <- as.matrix(all[, geno_cols])

alt_counts <- matrix(NA_integer_, nrow = nrow(geno_mat), ncol = ncol(geno_mat),
                     dimnames = dimnames(geno_mat))
gc()

alt_counts[geno_mat == "0/0"] <- 0
alt_counts[geno_mat == "0/1"] <- 1
alt_counts[geno_mat == "1/0"] <- 1
alt_counts[geno_mat == "1/1"] <- 2

sum_alt <- rowSums(alt_counts, na.rm = TRUE)
n_alleles <- rowSums(!is.na(alt_counts)) * 2
p <- (n_alleles - sum_alt) / n_alleles
q <- 1 - p

all$pi_FU <- 2 * p * q
all$Ho_FU <- rowSums(all[, geno_cols] == "0/1", na.rm = TRUE) / length(geno_cols)

piS <- mean(all$Ho_FU[all$cat1 == "SYN"], na.rm = TRUE)
PiS <- mean(all$pi_FU[all$cat1 == "SYN"], na.rm = TRUE)
piN <- mean(all$Ho_FU[all$cat1 == "MIS"], na.rm = TRUE)
PiN <- mean(all$pi_FU[all$cat1 == "MIS"], na.rm = TRUE)
piIG <- mean(all$Ho_FU[all$cat1 == "IG"], na.rm = TRUE)
PiIG <- mean(all$pi_FU[all$cat1 == "IG"], na.rm = TRUE)
piH <- mean(all$Ho_FU[all$cat1 == "LOF"], na.rm = TRUE)
PiH <- mean(all$pi_FU[all$cat1 == "LOF"], na.rm = TRUE)
PiA <- mean(all$pi_FU, na.rm = TRUE)
PiN_PiS <- PiN / PiS
PiH_Pis <- PiH / PiS

results_FU <- data.frame(
  pop = "FU",
  he_s = PiS,
  he_ns = PiN,
  he_h = PiH,
  he_ig = PiIG,
  he_all = PiA,
  ho_s = piS,
  ho_ns = piN,
  ho_h = piH,
  ho_ig = piIG
)

# Combine results
results <- rbind(results, results_WP, results_FU)

# Prepare data for plotting
df_He <- results %>%
  select(pop, starts_with("he_")) %>%
  pivot_longer(
    cols = starts_with("he_"),
    names_to = "category",
    values_to = "value"
  ) %>%
  mutate(category = recode(category,
                           "he_s" = "SYN",
                           "he_ns" = "MISSENSE",
                           "he_h" = "HIGH",
                           "he_ig" = "IG", 
                           "he_all" = "ALL"))

df_He$category[df_He$category == "MISSENSE"] <- "MIS"
df_He$category[df_He$category == "HIGH"] <- "LOF"
df_He$category <- factor(df_He$category, levels = c("ALL", "IG", "SYN", "MIS", "LOF"))
df_He$pop[df_He$pop == "SP"] <- "LSP"
df_He$pop[df_He$pop == "WP"] <- "LWP"
df_He$pop <- factor(df_He$pop, levels = c("LSP", "LWP", "FU"))
df_He$Category_label <- rep(c(81099, 53136, 764, 3943762, 1000000), 3)
#df_He$Category_label <- rep(c(2352,1772,22,84267, 1000000),3)

# Create line dataframe for Y-axis
population_names <- unique(df_He$pop)
first_population <- population_names[1]

line_df_median <- data.frame(
  x = rep(0, length(population_names)),
  xend = rep(0, length(population_names)),
  y = ifelse(population_names == first_population, -0.0002, NA),
  yend = ifelse(population_names == first_population, 0.2502, NA),
  pop = population_names
)

df_He <- df_He[df_He$category != "ALL", ]

# Plot expected heterozygosity
He_plot <- ggplot(df_He, aes(x = category, y = value, fill = category)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_segment(data = line_df_median, 
               aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE,
               color = "black", size = 0.5) +
  geom_text(aes(label = Category_label, y = value + 0.01),
            size = 6,
            vjust = 0, color = "black") +
  facet_wrap(~ pop) +
  scale_fill_manual(
    values = c("palevioletred1", "darkorange", "red", "darkred")
  ) +
  scale_y_continuous(
    limits = c(-0.01, 0.253),
    breaks = seq(0, 0.25, by = 0.05),
    labels = scales::label_number(accuracy = 0.01),
    expand = c(0, 0)
  ) +
  theme_minimal() +
  theme(
    strip.text = element_blank(),
    axis.title.y = element_text(size = 22, color = "black", margin = margin(r = 10)),
    axis.text.x = element_text(size = 22, color = "black"),
    axis.text.y = element_text(size = 22, color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks.y = element_line(color = "black", size = 0.8),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Expected Heterozygosity (He)",
    fill = NULL
  )

He_plot

save(df_He, He_plot, line_df_median, file = "He_plot.RData")
#save(df_He, He_plot, line_df_median, file = "Z_He_plot.RData")
save(piN_piS, pin_pis, file = "pin_pis_all.RData")
#save(piN_piS, pin_pis, file = "Z_pin_pis_all.RData")
