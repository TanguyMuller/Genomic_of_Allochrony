#   This script combines BayPass output files across multiple runs
#   to compute local scores (Lindley process) and visualize genome-wide
#   association signals.

library(data.table)
source("baypass_utils.R")

#-----------------------------------------------------------
# 1. Define input files and read contrast data
#-----------------------------------------------------------
fichiers <- sprintf("AsubGlRc_%d_summary_contrast.out.bz2", 1:100)
fichiers_pi <- sprintf("AsubGlRc_%d_summary_pi_xtx.out.bz2", 1:100)

Z <- fread("Z_gl_summary_contrast.out")
Z <- Z[Z$CONTRAST == "2", ]

#-----------------------------------------------------------
# 2. Function to read compressed BayPass output
#-----------------------------------------------------------
lire_fichier_bz2 <- function(fichier) {
  con <- bzfile(fichier, open = "r")
  df <- read.table(con, header = TRUE)
  close(con)
  return(df)
}

#-----------------------------------------------------------
# 3. Read and merge all summary files
#-----------------------------------------------------------
list_dfs <- lapply(fichiers, lire_fichier_bz2)
list_pi <- lapply(fichiers_pi, lire_fichier_bz2)

df_combine <- do.call(rbind, list_dfs)
df_combine <- df_combine[df_combine$CONTRAST == "2", ]

colnames(df_combine)[colnames(df_combine) == "log10.1.pval."] <- "log"
colnames(Z)[colnames(Z) == "log10(1/pval)"] <- "log"
df_combine <- rbind(df_combine, Z)

df_pi <- do.call(rbind, list_pi)
colnames(df_pi)[colnames(df_pi) == "log10.1.pval."] <- "log"

Z_pi <- fread("../baypass/anaZgl_all_summary_pi_xtx.out")
colnames(Z_pi)[colnames(Z_pi) == "log10(1/pval)"] <- "log"
df_pi <- rbind(df_pi, Z_pi)

#-----------------------------------------------------------
# 4. Add SNP positional information
#-----------------------------------------------------------
positions <- read.table("A.baypass.snpdet")
positions_Z <- fread("../baypass/Z.baypass.snpdet.bz2")
positions <- rbind(positions, positions_Z)

df_combine[, 1:2] <- positions[, 1:2]
df_combine$CONTRAST <- gsub("chr", "", df_combine$CONTRAST)
df_combine$CONTRAST <- factor(df_combine$CONTRAST, levels = c(as.character(1:49), "Z"), ordered = TRUE)

Z <- df_combine[df_combine$CONTRAST == "Z", ]

#-----------------------------------------------------------
# 5. Save results
#-----------------------------------------------------------
save.image(file = "baypass_combine.RData")
