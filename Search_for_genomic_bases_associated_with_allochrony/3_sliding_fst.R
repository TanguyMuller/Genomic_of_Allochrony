# Genomic scan analysis for chromosome Z
# Analysis of FST, nucleotide diversity (pi), Tajima's D, and read depth

library(poolfstat)
library(tidyr)
library(dplyr)
options(scipen = 999)

#-----------------------------------------------------------
# 1. Load and prepare pool-seq data
#-----------------------------------------------------------

# Pool names and sizes
pool.list <- as.character(c("01_LSP", "02_LWP", "03_FU", "04_VI", "05_VA",
                           "06_CA", "07_TA", "08_GR", "09_F1", "10_LateLSP"))

# Calculate haploid size for Z chromosome (sex chromosome adjustment)
diploide_size_VI <- 20
diploide_size_VA <- 26
diploide_size_F1 <- 17
m1 <- 0.32 / 0.5  # Male proportion
m2 <- 0.21 / 0.5
m3 <- 0.26 / 0.5

# Haploid size calculation (males have 1 Z, females have 2)
nbrVI <- m1 * diploide_size_VI * 2 + (1 - m1) * diploide_size_VI * 1
nbrVI <- round(nbrVI)
nbrVA <- m2 * diploide_size_VA * 2 + (1 - m2) * diploide_size_VA * 1
nbrVA <- round(nbrVA)
nbrF1 <- m3 * diploide_size_F1 * 2 + (1 - m3) * diploide_size_F1 * 1
nbrF1 <- round(nbrF1)

haploid_size_Z <- as.numeric(c(100, 100, 80, nbrVI, nbrVA, 70, 74, 70, nbrF1, 32))

# Load VCF data
data_Z <- vcf2pooldata("path/to/Z.pool.portugal.vcf.gz",
                       poolsizes = haploid_size_Z,
                       poolnames = pool.list,
                       min.cov.per.pool = 10)

#-----------------------------------------------------------
# 2. Calculate FST sliding windows between paires of populations
#-----------------------------------------------------------

SPWP <- pooldata.subset(data_Z, pool.index = c(1,2),min.cov.per.pool = 10, min.maf = 0.05)
SPLSP <- pooldata.subset(data_Z, pool.index = c(1,10),min.cov.per.pool = 10, min.maf = 0.05)
WPLSP <- pooldata.subset(data_Z, pool.index = c(2,10),min.cov.per.pool = 10, min.maf = 0.05)
FULSP <- pooldata.subset(data_Z, pool.index = c(3,10),min.cov.per.pool = 10, min.maf = 0.05)
SPFU <- pooldata.subset(data_Z, pool.index = c(1,3),min.cov.per.pool = 10, min.maf = 0.05)
WPFU<- pooldata.subset(data_Z, pool.index = c(2,3),min.cov.per.pool = 10, min.maf = 0.05)

# Calculate FST with sliding windows
fst1 <- computeFST(SPWP,sliding.window.size = 10)
fst2 <- computeFST(SPLSP,sliding.window.size = 10)
fst3 <- computeFST(WPLSP,sliding.window.size = 10)
fst4 <- computeFST(SPFU,sliding.window.size = 10)
fst5 <- computeFST(WPFU,sliding.window.size = 10)

save.image("fst_sliding_window.RData")

#-----------------------------------------------------------
# 3. Add mutation annotations for Figure 6
#-----------------------------------------------------------

fst_df <- data.frame(
  pos = fst1$sliding.windows.fvalues$CumMidPos,
  fst = fst1$sliding.windows.fvalues$MultiLocusFst
)

# Calculate FST per SNP
fst_region <- computeFST(SPWP)
fst_region <- data.frame(
  pos = SPWP@snp.info$Position,
  fst = fst_region$snp.Fstats$Fst
)

# Load mutation categories (TOLERANT, DELETERIOUS, LOF)
load("path/to/DEL_TOL_LOF_all.RData") # File of chr/position and cat 

# Filter for regions of interest
df_cat <- DEL_TOL_LOF[DEL_TOL_LOF$position > 12000000 & DEL_TOL_LOF$position < 14600000 | 
           DEL_TOL_LOF$position > 16600000 & DEL_TOL_LOF$position < 18800000, ]
# MIS variants were cut in TOLERANT and DELETERIOUS variants (moderate impact with no change in functional constraints and with a change in functional constraints, respectively)
df_cat <- df_cat[df_cat$cat == "TOLERANT" | 
                 df_cat$cat == "DELETERIOUS" | 
                 df_cat$cat == "LOF", ]

df_cat <- df_cat[, c(2, 58)]
colnames(df_cat) <- c("pos", "cat")
df_cat$pos <- as.numeric(df_cat$pos)

# Merge FST with mutation categories
fst_mis_region <- merge(x = df_cat, y = fst_region, by = "pos")
fst_mis_region$color <- ifelse(fst_mis_region$cat == "TOLERANT", "green", 
                               ifelse(fst_mis_region$cat == "DELETERIOUS", "red", "purple"))

# Sort by position
fst_mis_region <- fst_mis_region %>% arrange(pos)

# Separate by category
df_cat_fst_red <- fst_mis_region[fst_mis_region$color == "red", ]
df_cat_fst_green <- fst_mis_region[fst_mis_region$color == "green", ]
df_cat_fst_purple <- fst_mis_region[fst_mis_region$color == "purple", ]

#-----------------------------------------------------------
# 4. Calculate rolling mean for FST
#-----------------------------------------------------------

rolling_mean <- rollmean(zoo(fst_df$fst, order.by = fst_df$pos), 
                        k = 25, na.rm = TRUE)

rolling_df <- data.frame(
  pos = index(rolling_mean),
  mean = coredata(rolling_mean)
)

save.image("fst_with_mutations_sliding_fst.RData)
