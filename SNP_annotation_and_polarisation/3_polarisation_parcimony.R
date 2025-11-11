###############################
######## Polarization #########
###############################

# Load outgroup data
# Run 1_create_mpileup.sh and 2_launcher_mpileup2alleles.sh before
outgroup_compil <- read.delim("mpileup/load_ss_ES_outgroup.txt")

###################################
## Histogram of the polymorphism ##
###################################

# Function to process each group (OK, BO, WI)
process_outgroup <- function(prefix) {
  print(paste("Processing outgroup:", prefix))  # Ajout pour débogage
  
  allele_col <- paste0(prefix, "_allele")
  statut_col <- paste0(prefix, "_STATUT")
  c_col <- paste0(prefix, "_c")
  
  # Keep only polymorphic sites (Mult)
  polym <- outgroup_compil[outgroup_compil[[statut_col]] == "Mult", ]
  
  # Extract depths from allele field
  extract_values <- function(x) {
    if(is.na(x) || x == "") return(c(NA, NA))
    parts <- strsplit(x, ";")[[1]]
    if (length(parts) == 2) {
      values <- sapply(parts, function(part) {
        split_part <- strsplit(part, ":")[[1]]
        if(length(split_part) >= 2) {
          as.numeric(split_part[2])
        } else {
          NA
        }
      })
      return(values)
    } else {
      return(c(NA, NA))
    }
  }
  
  polym$extracted_values <- lapply(polym[[allele_col]], extract_values)
  polym$allele_1 <- sapply(polym$extracted_values, function(x) if (length(x) == 2) x[1] else NA)
  polym$allele_2 <- sapply(polym$extracted_values, function(x) if (length(x) == 2) x[2] else NA)
  
  # Filter: remove sites with total depth < 5
  polym$somme_alleles <- polym$allele_1 + polym$allele_2
  polym$allele_1[polym$somme_alleles < 5] <- NA
  polym$allele_2[polym$somme_alleles < 5] <- NA
  polym <- polym[!is.na(polym$allele_1) & !is.na(polym$allele_2), ]
  
  # Compute minor allele frequency
  polym$minor_freq <- pmin(polym$allele_1, polym$allele_2) / polym$somme_alleles
  
  # Compute 10% quantile threshold
  pdf(paste0("Hist_",prefix,".pdf"))
  hist(polym$minor_freq, breaks=200)
  dev.off()
  q10 <- quantile(polym$minor_freq, probs = 0.1, na.rm = TRUE)
  polym <- polym[polym$minor_freq < q10, ]
  
  # Extract base letters (e.g., A, T, C, G)
  get_bases <- function(x) {
    if(is.na(x) || x == "") return(c(NA, NA))
    parts <- strsplit(x, ";")[[1]]
    if (length(parts) == 2) {
      return(sapply(parts, function(part) {
        split_part <- strsplit(part, ":")[[1]]
        if(length(split_part) >= 1) {
          split_part[1]
        } else {
          NA
        }
      }))
    } else {
      return(c(NA, NA))
    }
  }
  
  bases <- t(sapply(polym[[allele_col]], get_bases))
  polym$allele_1_base <- bases[, 1]
  polym$allele_2_base <- bases[, 2]
  
  # Identify major allele (the one with highest depth)
  polym[[allele_col]] <- ifelse(polym$allele_1 >= polym$allele_2, polym$allele_1_base, polym$allele_2_base)
  polym[[statut_col]] <- "Uniq"
  
  # Update the main data frame
  idx <- outgroup_compil$SNPbi %in% polym$SNPbi
  outgroup_compil[idx, allele_col] <- polym[[allele_col]]
  outgroup_compil[idx, statut_col] <- polym[[statut_col]]
  
  print(paste("Processed", nrow(polym), "sites for", prefix))
  return(outgroup_compil)
}

# Apply the function to all three groups 
for (prefix in c("OK", "WI", "BO")) {
  result <- process_outgroup(prefix)
  if(!is.null(result)) {
    outgroup_compil <- result
  }
}
rm(result)

###########################################
###### Polarization from 3 outgroups ######
###########################################

# Step 1: filter DP < 5
preprocess_coverage <- function(df) {
  
  # Clone df
  result_df <- df
  
  # Applied if *_c < 5, then *_allele = NA and *_STATUT = "COVnul"
  
  # For OK
  low_cov_OK <- result_df$OK_c < 5 | is.na(result_df$OK_c)
  result_df$OK_allele[low_cov_OK] <- NA
  result_df$OK_STATUT[low_cov_OK] <- "COVnul"
  
  # For WI
  low_cov_WI <- result_df$WI_c < 5 | is.na(result_df$WI_c)
  result_df$WI_allele[low_cov_WI] <- NA
  result_df$WI_STATUT[low_cov_WI] <- "COVnul"
  
  # For BO
  low_cov_BO <- result_df$BO_c < 5 | is.na(result_df$BO_c)
  result_df$BO_allele[low_cov_BO] <- NA
  result_df$BO_STATUT[low_cov_BO] <- "COVnul"
  
  # Resume modifications
  n_ok_modified <- sum(low_cov_OK)
  n_bo_modified <- sum(low_cov_BO)
  n_wi_modified <- sum(low_cov_WI)
  
  cat(paste0("\nFilter on low cover sites (< 5):\n",
             "Modified sites for OK: ", n_ok_modified, " (", round(n_ok_modified/nrow(result_df)*100, 2), "%)\n",
             "Modified sites for BO: ", n_bo_modified, " (", round(n_bo_modified/nrow(result_df)*100, 2), "%)\n", 
             "Modified sites for WI: ", n_wi_modified, " (", round(n_wi_modified/nrow(result_df)*100, 2), "%)\n"))
  
  return(result_df)
}

preprocess_coverage(outgroup_compil)

# Step 2: Fonction of polarisation 
polarize_outgroup_fast <- function(outgroup_df) {
  
  # Clone df
  result_df <- outgroup_df
  
  # Condition 1: At least two statuts are "Uniq"
  result_df$nb_uniq <- rowSums(result_df[, c("OK_STATUT", "BO_STATUT", "WI_STATUT")] == "Uniq", na.rm = TRUE)
  result_df <- result_df[result_df$nb_uniq >= 2, ]
  
  OK_uniq <- result_df$OK_STATUT == "Uniq"
  BO_uniq <- result_df$BO_STATUT == "Uniq"
  WI_uniq <- result_df$WI_STATUT == "Uniq"
  
  # Initiate vector for ancestral allele
  result_df$ancestral_allele <- NA
  
  # Comparisons by paires
  # Case 1: OK and BO are both "Uniq" and they have the same allele 
  mask_OK_BO <- OK_uniq & BO_uniq & (result_df$OK_allele == result_df$BO_allele)
  result_df$ancestral_allele[mask_OK_BO] <- result_df$OK_allele[mask_OK_BO]
  
  # Case 2: OK and WI are both "Uniq" and they have the same allele
  mask_OK_WI <- OK_uniq & WI_uniq & (result_df$OK_allele == result_df$WI_allele) & is.na(result_df$ancestral_allele)
  result_df$ancestral_allele[mask_OK_WI] <- result_df$OK_allele[mask_OK_WI]
  
  # Case 3: BO and  WI are both "Uniq" and they have the same allele
  mask_BO_WI <- BO_uniq & WI_uniq & (result_df$BO_allele == result_df$WI_allele) & is.na(result_df$ancestral_allele)
  result_df$ancestral_allele[mask_BO_WI] <- result_df$BO_allele[mask_BO_WI]
  
  
  # Filtrer les lignes où un allèle ancestral a été déterminé
  result_df <- result_df[!is.na(result_df$ancestral_allele), ]
  
  # Comparer avec REF pour déterminer si REF est l'allèle ancestral
  result_df$REF_is_ANC <- ifelse(result_df$ancestral_allele == result_df$REF, 1, 0)
  
  # Calculer des statistiques
  nb_sites <- nrow(result_df)
  nb_ref_is_anc <- sum(result_df$REF_is_ANC)
  pourcentage <- round(nb_ref_is_anc / nb_sites * 100, 2)
  
  cat(paste0("\nRésultats de la polarisation:\n",
             "Nombre total de sites polarisés: ", nb_sites, "\n",
             "Sites où REF est l'allèle ancestral: ", nb_ref_is_anc, " (", pourcentage, "%)\n",
             "Sites où REF n'est PAS l'allèle ancestral: ", nb_sites - nb_ref_is_anc, 
             " (", round(100 - pourcentage, 2), "%)\n"))
  
  return(result_df)
}

# Fonction complète qui combine les deux étapes
analyze_polarization <- function(df) {
  cat("Étape 1: Prétraitement du filtre de couverture\n")
  preprocessed_df <- preprocess_coverage(df)
  
  cat("\nÉtape 2: Analyse de polarisation\n")
  polarized_df <- polarize_outgroup_fast(preprocessed_df)
  
  return(polarized_df)
}

# Exécution complète avec mesure du temps
system.time({
  result <- analyze_polarization(outgroup_compil)
})

# Sauvegarde du résultat
if(!is.null(result)) {
  write.table(result, file = "polarized_result.txt", 
              sep = "\t", quote = FALSE, row.names = FALSE)
}
