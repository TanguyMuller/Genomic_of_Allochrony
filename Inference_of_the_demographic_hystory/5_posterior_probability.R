# Script for ABC Random Forest parameter estimation with posterior distribution analysis

library(abcrf)
library(weights)
library(KScorrect)
library(writexl)
library(openxlsx)

# Load real data
reel <- read.table("real_data/obs_table_ref", h = T)

args <- commandArgs(trailingOnly = TRUE)
k <- as.integer(args[1])

name.list <- c("sim.Ne_ancpp", "sim.Ne_sp", "sim.Ne_wp", "sim.Ne_fu", "sim.Ne_sp_found", "sim.Ne_wp_found", "sim.Ne_fu_ancfound",
               "sim.Ne_wp_ancfound", "sim.Ne_fu_anc", "sim.Ne_wp_anc", "sim.Ne_bot_sp", "sim.Ne_bot_wp", "sim.Ne_bot_fu",
               "sim.tsplit_PP", "sim.tsplit_WP", "sim.t_bot", "sim.m_ancpp", "sim.m_spwp_anc", "sim.m_spwp_rec", "sim.m_spfu_anc",
               "sim.m_spfu_rec", "sim.m_wpfu_anc", "sim.m_wpfu_rec")

pdf(paste0("Posterior.scen_SP_WP_", name.list[k], ".pdf"), width = 300/72, height = 300/72)

predict.parameter <- data.frame("estimate_mean" = numeric(), "estimate_median" = numeric(), "90_CI_2.5" = numeric(), "90_CI_97.5" = numeric(),
                                "prior_NMAE_mean" = numeric(), "prior_NMEA_median" = numeric(),
                                "posterior_NMAE_mean" = numeric(), "posterior_NMEA_median" = numeric(),
                                "95_cov" = numeric())

df1 <- read.table("reftable_scen/reftable.recent_div.txt", h = T)

# Diagnostics: Display dimensions before selection
cat("Dimensions of real data:", dim(reel), "\n")
cat("Dimensions of reference table:", dim(df1), "\n")

reel_new <- reel[, c(1:32, 34, 36, 38:42, 44, 46, 48:52, 54, 56, 58:167)]
ref_table_1 <- as.data.frame(df1[, c(name.list[k], names(df1)[24:190])])
ref_table_1 <- ref_table_1[, c(1:33, 35, 37, 39:43, 45, 47, 49:53, 55, 57, 59:168)]
ref_table_1 <- na.omit(ref_table_1)

# Diagnostics: Check dimensions after selection
cat("Dimensions of reel_new:", dim(reel_new), "\n")
cat("Dimensions of ref_table_1:", dim(ref_table_1), "\n")
cat("Number of rows in ref_table_1 after na.omit:", nrow(ref_table_1), "\n")

# Diagnostics: Compare column names
cat("\n=== DIAGNOSTICS ===\n")
cat("Column names in ref_table_1 (excluding first column - response variable):\n")
predictors_train <- names(ref_table_1)[-1]
print(predictors_train)

cat("\nColumn names in reel_new:\n")
predictors_obs <- names(reel_new)
print(predictors_obs)

# Find missing columns
missing_in_obs <- setdiff(predictors_train, predictors_obs)
missing_in_train <- setdiff(predictors_obs, predictors_train)

if (length(missing_in_obs) > 0) {
    cat("\n!!! MISSING COLUMNS in reel_new:\n")
    print(missing_in_obs)
}

if (length(missing_in_train) > 0) {
    cat("\n!!! EXTRA COLUMNS in reel_new (not in ref_table_1):\n")
    print(missing_in_train)
}

# Correction: Ensure reel_new has exactly the same columns as ref_table_1 (except first column)
common_cols <- intersect(predictors_train, predictors_obs)
cat("\nNumber of common columns:", length(common_cols), "\n")

if (length(common_cols) < length(predictors_train)) {
    cat("\n!!! WARNING: Missing", length(predictors_train) - length(common_cols), "columns!\n")
    stop("Cannot continue: missing columns in observed data")
}

# Reorder reel_new to match ref_table_1 exactly
reel_new <- reel_new[, predictors_train]

cat("\nAfter correction - Dimensions of reel_new:", dim(reel_new), "\n")
cat("Verification: identical names?", identical(names(reel_new), predictors_train), "\n")

# Create regression model
formula <- paste("log10(", name.list[k], ") ~ .", sep = "")
formula <- as.formula(formula)
RFmodel <- regAbcrf(formula = formula,
                    data = ref_table_1,
                    ntree = 1000,
                    paral = TRUE)

# Plot 1: Prior vs posterior distribution
hist(log10(ref_table_1[, 1]),
     breaks = 100,
     col = "grey",
     ylim = c(0, 1.5),
     freq = FALSE,
     main = name.list[k],
     xlab = "",
     ylab = "probability density", cex.main = 1, cex.axis = 0.8)

posterior_RF <- predict(object = RFmodel,
                        obs = reel_new,
                        quantiles = c(0.05, 0.95),
                        training = ref_table_1,
                        paral = TRUE,
                        ncores = 32,
                        ntree = 1000,
                        rf.weights = TRUE,
                        post.err.med = TRUE)

wtd.hist(log10(ref_table_1[, 1]),
         breaks = 50,
         col = rgb(1, 1, 0, alpha = 0.5), freq = FALSE, add = TRUE,
         weight = posterior_RF$weights)

newdata <- data.frame("estimate_mean" = posterior_RF$expectation,
                      "estimate_median" = posterior_RF$med[1],
                      "90_CI_2.5" = posterior_RF$quantile[1],
                      "90_CI_97.5" = posterior_RF$quantile[2],
                      "prior_NMAE_mean" = posterior_RF$prior.NMAE.mean,
                      "prior_NMEA_median" = posterior_RF$prior.NMAE.med,
                      "posterior_NMAE_mean" = posterior_RF$post.NMAE.mean,
                      "posterior_NMEA_median" = posterior_RF$post.NMAE.med,
                      "95_cov" = posterior_RF$prior.coverage)

predict.parameter <- rbind(predict.parameter, newdata)

dev.off()

# Plot 2: Estimated vs Actual (Out-of-Bag) in PNG format
png(paste0("OOB_predictions_scen_SP_WP_", name.list[k], ".png"), 
    width = 2000, height = 2000, res = 300)

# Extract OOB predictions
oob_predictions <- RFmodel$model.rf$predictions
actual_values <- log10(ref_table_1[, 1])

# Create plot
plot(actual_values, oob_predictions,
     xlab = "Actual parameter (log10)",
     ylab = "Estimated parameter (log10)",
     main = paste("OOB predictions -", name.list[k]),
     pch = 16,
     col = rgb(0, 0, 0, alpha = 0.1),
     cex = 0.3)

# Add identity line (y=x)
abline(0, 1, col = "red", lwd = 2, lty = 2)

# Add linear regression line
abline(lm(oob_predictions ~ actual_values), col = "blue", lwd = 2)

# Calculate R²
r_squared <- cor(actual_values, oob_predictions)^2
legend("topleft", 
       legend = c(paste("R² =", round(r_squared, 3)),
                  "Identity line",
                  "Linear fit"),
       col = c("black", "red", "blue"),
       lty = c(0, 2, 1),
       lwd = c(0, 2, 2),
       pch = c(16, NA, NA),
       bty = "n",
       cex = 1.2)

dev.off()

predict.parameter$parameter <- name.list[k]
chemin_fichier_excel <- paste0("parameter_scen_SP_WP_", name.list[k], ".xlsx")
write_xlsx(predict.parameter, chemin_fichier_excel)


