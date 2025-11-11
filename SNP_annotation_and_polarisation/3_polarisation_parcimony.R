#!/usr/bin/env Rscript

# Polarization analysis using three outgroups
# Run 1_create_mpileup.sh and 2_launcher_mpileup2alleles.sh before running this script.

outgroup_compil <- read.delim("mpileup/load_outgroup.txt")
#outgroup_compil <- read.delim("mpileup/Z_load_outgroup.txt")

# ----------------------------
# Function to process outgroups
# ----------------------------
process_outgroup <- function(prefix) {
  cat("Processing outgroup:", prefix, "\n")
  
  allele_col <- paste0(prefix, "_allele")
  statut_col <- paste0(prefix, "_STATUT")
  c_col <- paste0(prefix, "_c")
  
  # Keep only polymorphic sites
  polym <- outgroup_compil[outgroup_compil[[statut_col]] == "Mult", ]
  
  # Extract allele depths
  extract_values <- function(x) {
    if (is.na(x) || x == "") return(c(NA, NA))
    parts <- strsplit(x, ";")[[1]]
    if (length(parts) == 2) {
      sapply(parts, function(part) {
        vals <- strsplit(part, ":")[[1]]
        if (length(vals) >= 2) as.numeric(vals[2]) else NA
      })
    } else c(NA, NA)
  }
  
  polym$values <- lapply(polym[[allele_col]], extract_values)
  polym$allele_1 <- sapply(polym$values, `[`, 1)
  polym$allele_2 <- sapply(polym$values, `[`, 2)
  polym$total_depth <- polym$allele_1 + polym$allele_2
  
  # Filter sites with total depth < 5
  polym <- polym[polym$total_depth >= 5, ]
  
  # Compute minor allele frequency and keep lowest 10%
  polym$minor_freq <- pmin(polym$allele_1, polym$allele_2) / polym$total_depth
  pdf(paste0("Hist_", prefix, ".pdf"))
  hist(polym$minor_freq, breaks = 200)
  dev.off()
  q10 <- quantile(polym$minor_freq, 0.1, na.rm = TRUE)
  polym <- polym[polym$minor_freq < q10, ]
  
  # Extract allele bases
  get_bases <- function(x) {
    if (is.na(x) || x == "") return(c(NA, NA))
    parts <- strsplit(x, ";")[[1]]
    if (length(parts) == 2)
      sapply(parts, function(p) strsplit(p, ":")[[1]][1])
    else c(NA, NA)
  }
  
  bases <- t(sapply(polym[[allele_col]], get_bases))
  polym$allele_1_base <- bases[, 1]
  polym$allele_2_base <- bases[, 2]
  
  # Keep the major allele (higher depth)
  polym[[allele_col]] <- ifelse(polym$allele_1 >= polym$allele_2,
                                polym$allele_1_base, polym$allele_2_base)
  polym[[statut_col]] <- "Uniq"
  
  # Update main data frame
  idx <- outgroup_compil$SNPbi %in% polym$SNPbi
  outgroup_compil[idx, allele_col] <- polym[[allele_col]]
  outgroup_compil[idx, statut_col] <- polym[[statut_col]]
  
  cat("Processed", nrow(polym), "sites for", prefix, "\n")
  outgroup_compil
}

# Apply function to all three outgroups
for (prefix in c("OK", "WI", "BO")) {
  outgroup_compil <- process_outgroup(prefix)
}

# -------------------------------
# Coverage filtering (< 5 reads)
# -------------------------------
preprocess_coverage <- function(df) {
  result_df <- df
  
  for (prefix in c("OK", "WI", "BO")) {
    cov_col <- paste0(prefix, "_c")
    allele_col <- paste0(prefix, "_allele")
    statut_col <- paste0(prefix, "_STATUT")
    low_cov <- result_df[[cov_col]] < 5 | is.na(result_df[[cov_col]])
    result_df[[allele_col]][low_cov] <- NA
    result_df[[statut_col]][low_cov] <- "COVnul"
    cat("Filtered", sum(low_cov), "low-coverage sites for", prefix, "\n")
  }
  result_df
}

# ------------------------
# Polarization procedure
# ------------------------
polarize_outgroup <- function(df) {
  result_df <- df
  result_df$nb_uniq <- rowSums(result_df[, c("OK_STATUT", "BO_STATUT", "WI_STATUT")] == "Uniq", na.rm = TRUE)
  result_df <- result_df[result_df$nb_uniq >= 2, ]
  
  OK_uniq <- result_df$OK_STATUT == "Uniq"
  BO_uniq <- result_df$BO_STATUT == "Uniq"
  WI_uniq <- result_df$WI_STATUT == "Uniq"
  result_df$ancestral_allele <- NA
  
  mask_OK_BO <- OK_uniq & BO_uniq & (result_df$OK_allele == result_df$BO_allele)
  mask_OK_WI <- OK_uniq & WI_uniq & (result_df$OK_allele == result_df$WI_allele)
  mask_BO_WI <- BO_uniq & WI_uniq & (result_df$BO_allele == result_df$WI_allele)
  
  result_df$ancestral_allele[mask_OK_BO] <- result_df$OK_allele[mask_OK_BO]
  result_df$ancestral_allele[mask_OK_WI & is.na(result_df$ancestral_allele)] <- result_df$OK_allele[mask_OK_WI]
  result_df$ancestral_allele[mask_BO_WI & is.na(result_df$ancestral_allele)] <- result_df$BO_allele[mask_BO_WI]
  
  result_df <- result_df[!is.na(result_df$ancestral_allele), ]
  result_df$REF_is_ANC <- as.integer(result_df$ancestral_allele == result_df$REF)
  
  nb_sites <- nrow(result_df)
  nb_ref_is_anc <- sum(result_df$REF_is_ANC)
  cat("Polarization results:\n",
      "Total polarized sites:", nb_sites, "\n",
      "REF = ancestral:", nb_ref_is_anc, "(", round(nb_ref_is_anc / nb_sites * 100, 2), "%)\n",
      "REF ≠ ancestral:", nb_sites - nb_ref_is_anc, "(", round(100 - (nb_ref_is_anc / nb_sites * 100), 2), "%)\n")
  
  result_df
}

# ------------------------
# Full pipeline execution
# ------------------------
analyze_polarization <- function(df) {
  cat("Step 1: Coverage filtering\n")
  pre_df <- preprocess_coverage(df)
  cat("\nStep 2: Polarization\n")
  polarize_outgroup(pre_df)
}

system.time({
  result <- analyze_polarization(outgroup_compil)
})

if (!is.null(result)) {
  write.table(result, "polarized_result.txt", sep = "\t", quote = FALSE, row.names = FALSE)
}
