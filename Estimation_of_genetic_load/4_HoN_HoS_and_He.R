### Nouveau HeN/HeS
library(tidyverse)
library(dplyr)

# Load the GT filte
all <- read.delim("load.GT.txt", skip=1, header=FALSE)
all <- all[,c(1:43,81:94,105:107)]
SP <- paste0("SP_", 1:17)
WP <- paste0("WP_", 1:12)
FU <- paste0("FU_", 1:10)
SP_2 <- paste0("SP_", 18:25)
WP_2 <- paste0("WP_", 13:18)
colnames(all) <- c("Chromosome", 
                      "position", 
                      "REF", 
                      "ALT", 
                      SP, 
                      WP,
                      FU,
                      SP_2,
                      WP_2,
                      "SNPEFF_2", 
                      "SNPEFF_3", 
                      "SNPEFF_11")
all$SNPEFF_2 <- as.factor(all$SNPEFF_2)
all$SNPEFF_3 <- as.factor(all$SNPEFF_3)
all$SNPEFF_11 <- as.factor(all$SNPEFF_11)

# Add columns for region information
# Adding column with amino acid names when applicable:
all$AA1 <- as.factor(str_sub(all$SNPEFF_11, 3, 5)) # 
all$AA2 <- as.factor(ifelse(all$SNPEFF_3=="HIGH", NA, str_sub(all$SNPEFF_11, -3, -1)))

# Adding column with simplified "exon" / "not exon" status:
all$exon_position_status <- as.factor(ifelse(is.na(all$SNPEFF_11)==TRUE, "not_exon", "EXON"))

# Adding another column with simplified "intron" / "not intron" status:
all$intron_position_status <- as.factor(ifelse(all$SNPEFF_2=="intron_variant", "INTRON", "not_intron"))

# Adding another column with simplified "intergenic" / "not intergenic" status:
all$intergenic_position_status <- as.factor(ifelse(all$SNPEFF_2=="intergenic_region", "INTERGENIC", "not_intergenic"))

all <- all[!is.na(all$intergenic_position_status), ]

all$SNPbi <- paste(all$Chromosome,"__",all$position, sep="")

# -- Storage structure at each step i, number of kept SNPs, number of skipped SNPs
i<-0
v<-rep(0,2)
logs <- data.frame( filt_id=v, filt_name=v, remain=v, discard=v  )
rm(v)

# -- data 
logs[i+1,]<- c(i, "start", length(all$Chromosome) , 0 )
logs[i+1,]

# First I discard duplicate and triplicate site
i<-i+1
nb_discarded <- 0
nb_discarded <- length(all[duplicated(all$SNPbi),]$SNPbi)
all <- all[!duplicated(all$SNPbi),]
logs[i+1,]<- c(i, "duplicate and triplicate", length(all$SNPbi) , nb_discarded )
logs[i+1,]

# "N" genotypes filter
i<-i+1
nb_discarded <- 0
nb_discarded <- length(all[all$REF=="N",]$SNPbi )
all <- all[all$REF!="N",]
nb_discarded <- nb_discarded + length(all[all$ALT=="N",]$SNPbi )
all <- all[all$ALT!="N",]
logs[i+1,]<- c(i, "REF(/ALT)==N", length(all$SNPbi) , nb_discarded )
logs[i+1,]

# Amino Acid filter
i<-i+1
nb_discarded <- 0
# We discard SNPs with "???" instead of AA1, and "t1?" instead of AA2 :
nb_discarded <- length(na.omit(all[ all$AA1=="???", ]$SNPbi)) + length(na.omit(all[ all$AA2=="t1?", ]$SNPbi))
all <- subset(all,!(all$AA1 == "???") | is.na(all$AA1))
all <- subset(all,!(all$AA2 == "t1?") | is.na(all$AA2))
logs[i+1,]<- c(i, "AA not identified", length(all$SNPbi) , nb_discarded )
logs[i+1,]

## On ajoute une colonne avec les catégorie
all$cat1 <- "other" 
all[ !is.na(all$intergenic_position_status) & all$intergenic_position_status=="INTERGENIC" , ]$cat1 <- "IG"
all[ all$SNPEFF_3=="LOW" & all$exon_position_status=="EXON", ]$cat1 <- "SYN"
all[ all$SNPEFF_3=="MODERATE", ]$cat1 <- "MIS"
all[ all$SNPEFF_3=="HIGH" & all$exon_position_status=="EXON", ]$cat1 <- "LOF"
table(all$cat1)

#
### PiN/PiS
# Créer un dataframe vide pour stocker les résultats
results_SP <- data.frame()

# Itérer à travers chaque individu de SP1 à SP25
for (i in 1:25) {
  # Sélectionner la colonne correspondant à l'individu
  colname <- paste0("SP_", i)
  
  
  # Sélectionner les catégories
  NS <- all[!is.na(all$cat1) & all$cat1 == "MIS", ]
  Ho_NS <- sum(NS[[colname]] == "0/1") / nrow(NS)
  S <- all[!is.na(all$cat1) & all$cat1 == "SYN", ]
  Ho_S <- sum(S[[colname]] == "0/1") / nrow(S)
  NSH <- all[!is.na(all$cat1) & all$cat1 == "LOF", ]
  Ho_NSH <- sum(NSH[[colname]] == "0/1") /nrow(NSH)
  
  
  # Créer une ligne pour l'individu
  individual_results <- data.frame(
    ind = colname,
    sum_PiS_Ho = Ho_S,
    sum_PiN_Ho = Ho_NS,
    sum_PiNSH_Ho = Ho_NSH,
    PiN_PiS_Ho = ifelse(Ho_S != 0, Ho_NS / Ho_S, NA),
    PiNH_PiS_Ho = ifelse(Ho_S != 0, Ho_NSH / Ho_S, NA)
  )
  
  # Ajouter les résultats pour cet individu au dataframe
  results_SP <- rbind(results_SP, individual_results)
}

# Afficher les résultats
print(results_SP)


# Calcul des statistiques pour la population "WP"
results_WP <- data.frame()

# Itérer à travers chaque individu de WP1 à WP25
for (i in 1:18) {
  # Sélectionner la colonne correspondant à l'individu
  colname <- paste0("WP_", i)
  
  
  # Sélectionner les catégories
  NS <- all[!is.na(all$cat1) & all$cat1 == "MIS", ]
  Ho_NS <- sum(NS[[colname]] == "0/1") / nrow(NS)
  S <- all[!is.na(all$cat1) & all$cat1 == "SYN", ]
  Ho_S <- sum(S[[colname]] == "0/1") / nrow(S)
  NSH <- all[!is.na(all$cat1) & all$cat1 == "LOF", ]
  Ho_NSH <- sum(NSH[[colname]] == "0/1") /nrow(NSH)
  
  
  # Créer une ligne pour l'individu
  individual_results <- data.frame(
    ind = colname,
    sum_PiS_Ho = Ho_S,
    sum_PiN_Ho = Ho_NS,
    sum_PiNSH_Ho = Ho_NSH,
    PiN_PiS_Ho = ifelse(Ho_S != 0, Ho_NS / Ho_S, NA),
    PiNH_PiS_Ho = ifelse(Ho_S != 0, Ho_NSH / Ho_S, NA)
  )
  
  # Ajouter les résultats pour cet individu au dataframe
  results_WP <- rbind(results_WP, individual_results)
}

# Afficher les résultats
print(results_WP)


# Calcul des statistiques pour la population "FU"
results_FU <- data.frame()

# Itérer à travers chaque individu de FU1 à FU25
for (i in 1:10) {
  # Sélectionner la colonne correspondant à l'individu
  colname <- paste0("FU_", i)
  
  
  # Sélectionner les catégories
  NS <- all[!is.na(all$cat1) & all$cat1 == "MIS", ]
  Ho_NS <- sum(NS[[colname]] == "0/1") / nrow(NS)
  S <- all[!is.na(all$cat1) & all$cat1 == "SYN", ]
  Ho_S <- sum(S[[colname]] == "0/1") / nrow(S)
  NSH <- all[!is.na(all$cat1) & all$cat1 == "LOF", ]
  Ho_NSH <- sum(NSH[[colname]] == "0/1") /nrow(NSH)
  
  
  # Créer une ligne pour l'individu
  individual_results <- data.frame(
    ind = colname,
    sum_PiS_Ho = Ho_S,
    sum_PiN_Ho = Ho_NS,
    sum_PiNSH_Ho = Ho_NSH,
    PiN_PiS_Ho = ifelse(Ho_S != 0, Ho_NS / Ho_S, NA),
    PiNH_PiS_Ho = ifelse(Ho_S != 0, Ho_NSH / Ho_S, NA)
  )
  
  # Ajouter les résultats pour cet individu au dataframe
  results_FU <- rbind(results_FU, individual_results)
}

# Afficher les résultats
print(results_FU)

# Fusionner les résultats
piN_piS_Z <- rbind(results_SP,results_WP,results_FU)
piN_piS_Z$pop<-c(rep("LSP",25), rep("LWP",18), rep("FU",10))
piN_piS<-piN_piS_Z

piN_piS$pop <- factor(piN_piS$pop, levels = c("LSP", "LWP", "FU"))
color_pinpis<-c("brown1", "cadetblue3", "chartreuse3")

pin_pis <- ggplot(data=piN_piS, aes(pop, PiN_PiS_Ho, fill=pop))+
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.8), 
              size = 2, color = "black", alpha = 0.5) +
  scale_fill_manual(values = color_pinpis)+
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
  labs(
    x = NULL,
    y = "NS/S Ratio"
  )

pin_pis

save(pin_pis, piN_piS, file = "pin_pis_all.RData")

#test stat
load("pin_pis_all.RData")
shapiro_test_results <- by(piN_piS$PiN_PiS_Ho, interaction(piN_piS$pop), shapiro.test)
shapiro_test_results
ggplot(piN_piS, aes(sample = PiN_PiS_Ho)) +
  facet_wrap(~ pop ) +
  stat_qq() +
  stat_qq_line() +
  theme_minimal() +
  labs(title = "QQ Plots by Population and Chromosome Type",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")
kruskal.test(PiN_PiS_Ho ~ pop, data = piN_piS)
library(FSA)
dunnTest(PiN_PiS_Ho ~ pop, data = piN_piS, method = "bonferroni")

#### het Ho and He
# He
geno_cols <- grep("^SP_", colnames(all), value = TRUE)
geno_mat <- as.matrix(all[, geno_cols])

alt_counts <- matrix(NA_integer_, nrow = nrow(geno_mat), ncol = ncol(geno_mat),
                     dimnames = dimnames(geno_mat))
gc()
alt_counts[geno_mat == "0/0"] <- 0
alt_counts[geno_mat == "0/1"] <- 1
alt_counts[geno_mat == "1/0"] <- 1  # si présent
alt_counts[geno_mat == "1/1"] <- 2
# Somme des allèles alt par site
sum_alt <- rowSums(alt_counts, na.rm = TRUE)
# Nombre total d’allèles (2 * nb d’individus sans NA)
n_alleles <- rowSums(!is.na(alt_counts)) * 2
# Fréquence allèle alt
p <- (n_alleles-sum_alt) / n_alleles
q <- 1 - p
# Pi = 2pq
all$pi_SP <- 2 * p * q
all$Ho_SP <- rowSums(all[, geno_cols] == "0/1", na.rm = TRUE) / length(geno_cols)
piS <- mean(all$Ho_SP[all$cat1 == "SYN"], na.rm = TRUE)
PiS <- mean(all$pi_SP[all$cat1 == "SYN"], na.rm = TRUE)       # PiS
piN <- mean(all$Ho_SP[all$cat1 == "MIS"], na.rm = TRUE)
PiN <- mean(all$pi_SP[all$cat1 == "MIS"], na.rm = TRUE)  # PiN
piIG <- mean(all$Ho_SP[all$cat1 == "IG"], na.rm = TRUE)
PiIG <- mean(all$pi_SP[all$cat1 == "IG"], na.rm = TRUE)
piH <- mean(all$Ho_SP[all$cat1 == "LOF"], na.rm = TRUE)
PiH <- mean(all$pi_SP[all$cat1 == "LOF"], na.rm = TRUE)
PiA <- mean(all$pi_SP, na.rm = TRUE)
PiN_PiS <- PiN / PiS
PiH_Pis <- PiH / PiS

hist(all$pi_SP[all$cat1 == "SYN"], na.rm = TRUE, breaks = 20)
hist(all$Ho_SP[all$cat1 == "SYN"], na.rm = TRUE, breaks = 100)
hist(all$pi_SP[all$cat1 == "MIS"], na.rm = TRUE, breaks = 20)
hist(all$Ho_SP[all$cat1 == "MIS"], na.rm = TRUE, breaks = 100)

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

# WP
geno_cols <- grep("^WP_", colnames(all), value = TRUE)
geno_mat <- as.matrix(all[, geno_cols])

alt_counts <- matrix(NA_integer_, nrow = nrow(geno_mat), ncol = ncol(geno_mat),
                     dimnames = dimnames(geno_mat))
gc()
alt_counts[geno_mat == "0/0"] <- 0
alt_counts[geno_mat == "0/1"] <- 1
alt_counts[geno_mat == "1/0"] <- 1  # si présent
alt_counts[geno_mat == "1/1"] <- 2
# Somme des allèles alt par site
sum_alt <- rowSums(alt_counts, na.rm = TRUE)
# Nombre total d’allèles (2 * nb d’individus sans NA)
n_alleles <- rowSums(!is.na(alt_counts)) * 2
# Fréquence allèle alt
p <- (n_alleles-sum_alt) / n_alleles
q <- 1 - p
# Pi = 2pq
all$pi_WP <- 2 * p * q
all$Ho_WP <- rowSums(all[, geno_cols] == "0/1", na.rm = TRUE) / length(geno_cols)
piS <- mean(all$Ho_WP[all$cat1 == "SYN"], na.rm = TRUE)
PiS <- mean(all$pi_WP[all$cat1 == "SYN"], na.rm = TRUE)       # PiS
piN <- mean(all$Ho_WP[all$cat1 == "MIS"], na.rm = TRUE)
PiN <- mean(all$pi_WP[all$cat1 == "MIS"], na.rm = TRUE)  # PiN
piIG <- mean(all$Ho_WP[all$cat1 == "IG"], na.rm = TRUE)
PiIG <- mean(all$pi_WP[all$cat1 == "IG"], na.rm = TRUE)
piH <- mean(all$Ho_WP[all$cat1 == "LOF"], na.rm = TRUE)
PiH <- mean(all$pi_WP[all$cat1 == "LOF"], na.rm = TRUE)
PiA <- mean(all$pi_WP, na.rm = TRUE)
PiN_PiS <- PiN / PiS
PiH_Pis <- PiH / PiS

hist(all$pi_WP[all$cat1 == "SYN"], na.rm = TRUE, breaks = 20)
hist(all$Ho_WP[all$cat1 == "SYN"], na.rm = TRUE, breaks = 100)
hist(all$pi_WP[all$cat1 == "MIS"], na.rm = TRUE, breaks = 20)
hist(all$Ho_WP[all$cat1 == "MIS"], na.rm = TRUE, breaks = 100)

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

#FU
geno_cols <- grep("^FU_", colnames(all), value = TRUE)
geno_mat <- as.matrix(all[, geno_cols])

alt_counts <- matrix(NA_integer_, nrow = nrow(geno_mat), ncol = ncol(geno_mat),
                     dimnames = dimnames(geno_mat))
gc()
alt_counts[geno_mat == "0/0"] <- 0
alt_counts[geno_mat == "0/1"] <- 1
alt_counts[geno_mat == "1/0"] <- 1  # si présent
alt_counts[geno_mat == "1/1"] <- 2
# Somme des allèles alt par site
sum_alt <- rowSums(alt_counts, na.rm = TRUE)
# Nombre total d’allèles (2 * nb d’individus sans NA)
n_alleles <- rowSums(!is.na(alt_counts)) * 2
# Fréquence allèle alt
p <- (n_alleles-sum_alt) / n_alleles
q <- 1 - p
# Pi = 2pq
all$pi_FU <- 2 * p * q
all$Ho_FU <- rowSums(all[, geno_cols] == "0/1", na.rm = TRUE) / length(geno_cols)
piS <- mean(all$Ho_FU[all$cat1 == "SYN"], na.rm = TRUE)
PiS <- mean(all$pi_FU[all$cat1 == "SYN"], na.rm = TRUE)       # PiS
piN <- mean(all$Ho_FU[all$cat1 == "MIS"], na.rm = TRUE)
PiN <- mean(all$pi_FU[all$cat1 == "MIS"], na.rm = TRUE)  # PiN
piIG <- mean(all$Ho_FU[all$cat1 == "IG"], na.rm = TRUE)
PiIG <- mean(all$pi_FU[all$cat1 == "IG"], na.rm = TRUE)
piH <- mean(all$Ho_FU[all$cat1 == "LOF"], na.rm = TRUE)
PiH <- mean(all$pi_FU[all$cat1 == "LOF"], na.rm = TRUE)
PiA <- mean(all$pi_FU, na.rm = TRUE)
PiN_PiS <- PiN / PiS
PiH_Pis <- PiH / PiS

hist(all$pi_FU[all$cat1 == "SYN"], na.rm = TRUE, breaks = 20)
hist(all$Ho_FU[all$cat1 == "SYN"], na.rm = TRUE, breaks = 100)
hist(all$pi_FU[all$cat1 == "MIS"], na.rm = TRUE, breaks = 20)
hist(all$Ho_FU[all$cat1 == "MIS"], na.rm = TRUE, breaks = 100)

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
results <- rbind(results, results_WP, results_FU)

# Figure 
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
df_He$Category_label <- rep(c(81099,53136,764,3943762, 1000000),3)

# Créer un dataframe pour la ligne Y uniquement pour la première Population
# Récupérer toutes les populations
population_names <- unique(df_He$pop)
first_population <- population_names[1]
line_df_median <- data.frame(
  x = rep(0, length(population_names)),
  xend = rep(0, length(population_names)),
  y = ifelse(population_names == first_population, -0.0002, NA),
  yend = ifelse(population_names == first_population, 0.2502, NA),
  pop = population_names  # Assurez-vous que le nom de la colonne correspond à celui utilisé dans df_He
)
df_He <- df_He[df_He$category!="ALL",]
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
    breaks = seq(0, 0.25, by = 0.05),       # 👉 Ticks tous les 0.05
    labels = scales::label_number(accuracy = 0.01),  # 👉 Affichage à deux décimales
    expand = c(0, 0)
  )+
  theme_minimal() +
  theme(
    strip.text = element_blank(),
    axis.title.y = element_text(size = 22, color = "black", margin = margin(r = 10)),  # Titre axe Y plus grand
    axis.text.x = element_text(size = 22, color = "black"),
    axis.text.y = element_text(size = 22, color = "black"),  # Valeurs Y plus grosses
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

save(df_He,He_plot,line_df_median,file = "He_plot.RData")
save(piN_piS,pin_pis, file = "pin_pis_all.RData")

#test stats
all_sp <- all[,c(67,68)]
all_wp <- all[,c(67,70)]
all_fu <- all[,c(67,72)]

shapiro_test_results <- by(all_sp$pi_SP, interaction(all_sp$cat1), shapiro.test)
kruskal.test(pi_SP ~ cat1, data = all_sp)
library(FSA)
dunnTest(pi_SP ~ cat1, data = all_sp, method = "bonferroni")
