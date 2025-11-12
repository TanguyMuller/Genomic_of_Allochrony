#   This script read multiple gzipped mosdepth coverage files,
#   combines them by group (LSP / LWP), and computes
#   median coverage ratios per window for the Figure 6.

# Set working directory
setwd("/path/mosdepth/")
library(dplyr)
library(purrr)

# List all .gz files
fichiers <- list.files(pattern = "gz", full.names = TRUE)
n_fichiers <- length(fichiers)

# Read each file and store in separate data frames (df_1, df_2, ...)
for (i in 1:n_fichiers) {
  df <- read.table(gzfile(fichiers[i]), header = FALSE, sep = "\t")
  assign(paste0("df_", i), df)
}

# ------------------------------------------------
# Function to combine and summarize mosdepth tables
# ------------------------------------------------
combine_dfs <- function(df_list) {
  # Add source ID to each dataframe
  df_list_with_id <- purrr::map(1:length(df_list), ~dplyr::mutate(df_list[[.x]], source_id = .x))
  
  # Combine all into one dataframe
  combined_df <- dplyr::bind_rows(df_list_with_id)
  
  # Compute median of V4 (coverage) per window
  result <- combined_df %>%
    dplyr::group_by(V2, V3) %>%
    dplyr::summarize(V4_mean = median(V4, na.rm = TRUE), .groups = "drop")
  
  return(result)
}

# ------------------------------------------------
# Combine SP and WP samples separately
# ------------------------------------------------

# SP group
df_names_SP <- paste0("df_", c(1:25))
df_list_SP <- mget(df_names_SP, envir = .GlobalEnv)
result_df_SP <- combine_dfs(df_list_SP)

# WP group
df_names_WP <- paste0("df_", c(26:44))
df_list_WP <- mget(df_names_WP, envir = .GlobalEnv)
result_df_WP <- combine_dfs(df_list_WP)

# ------------------------------------------------
# Compute median coverages and ratio between groups
# ------------------------------------------------
SP_med <- median(result_df_SP$V4_mean)
WP_med <- median(result_df_WP$V4_mean)

# Ratio of normalized coverage difference per window
result_df_SP$ratio <- (result_df_SP$V4_mean - result_df_WP$V4_mean) / ((SP_med + WP_med) / 2)
