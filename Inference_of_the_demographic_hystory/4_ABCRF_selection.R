# Script for ABC-RF Model Selection, Parameter Estimation, and Visualization

library(abcrf)
library(weights)
library(KScorrect)
library(writexl)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyr)

# Load real data
real_data <- read.table("real_data/obs_table_ref", h=T)

df1 <- read.table("reftable_scen/reftable.concomitent_div.txt", h=T)
df1 <- na.omit(df1)
df4 <- read.table("reftable_scen/reftable.recent_div.txt", h=T)
df4 <- na.omit(df4)
df7 <- read.table("reftable_scen/reftable.old_div.txt", h=T)
df7 <- na.omit(df7)

# Create folders for PCA and OOB plots
if (!dir.exists("PCA_plots")) {
  dir.create("PCA_plots")
}
if (!dir.exists("OOB_plots")) {
  dir.create("OOB_plots")
}

# Initialize lists to store results
confusion_matrices <- list()
prior_errors <- list()
model_selection_results_votes1 <- list()
model_selection_results_votes2 <- list()
model_selection_results_votes3 <- list()
model_selection_results_posterior <- list()

# Number of simulations and repeats
num_of_sims = 10000
num_repeats = 10

# Final dataframe to store results
df_final <- data.frame(
  subset_data = integer(0),
  global_error_rate = numeric(0),
  sd_prior_error = numeric(0),
  mean_votes1 = numeric(0),
  mean_votes2 = numeric(0),
  mean_votes3 = numeric(0),
  posterior_probability = numeric(0),
  local_error_rate = numeric(0),
  sd_votes1 = numeric(0),
  sd_votes2 = numeric(0),
  sd_votes3 = numeric(0),
  sd_posterior_probability = numeric(0),
  selected_model = character(0)
)

# Loop to repeat the operation 10 times
for (j in 1:num_repeats) {
  
  # Random subsampling
  set.seed(j)
  df1_sampled <- df1[sample(nrow(df1), num_of_sims), ]
  df4_sampled <- df4[sample(nrow(df4), num_of_sims), ]
  df7_sampled <- df7[sample(nrow(df7), num_of_sims), ]
  
  # Create the model factor
  model <- as.factor(c(rep("concomitent_div", num_of_sims),
                       rep("recent_div", num_of_sims),
                       rep("old_div", num_of_sims)))
  
  # Summary statistics
  sumstats <- rbind(df1_sampled[, 19:185], df4_sampled[, 24:190], df7_sampled[, 24:190])
  # We remove some summary statistics that didn't work in our case 
  sumstats <- sumstats[, c(1:32, 34, 36, 38:42, 44, 46, 48:52, 54, 56, 58:167)]
  sumstats$random <- rnorm(n = nrow(sumstats), mean = 0, sd = 1)
  real_new <- real_data[, c(1:32, 34, 36, 38:42, 44, 46, 48:52, 54, 56, 58:167)]
  real_new$random <- rnorm(n = nrow(real_new), mean = 0, sd = 1)
  
  # Create the reference table
  ref_table <- cbind(model, sumstats)
  ref_table <- na.omit(ref_table)
  
  # ========================================================================
  # PCA Plots
  # ========================================================================
  
  scen1_data <- sumstats[model == "concomitent_div", ]
  scen2_data <- sumstats[model == "recent_div", ]
  scen3_data <- sumstats[model == "old_div", ]
  
  pca1 <- prcomp(scen1_data, scale. = TRUE, center = TRUE)
  obs_proj1 <- predict(pca1, newdata = real_new)
  
  pca1_df <- data.frame(PC1 = pca1$x[, 1], PC2 = pca1$x[, 2], Type = "Simulated")
  obs1_df <- data.frame(PC1 = obs_proj1[, 1], PC2 = obs_proj1[, 2], Type = "Observed")
  plot1_data <- rbind(pca1_df, obs1_df)
  
  p1 <- ggplot(plot1_data, aes(x = PC1, y = PC2)) +
    geom_point(data = subset(plot1_data, Type == "Simulated"), 
               aes(color = Type), alpha = 0.6, size = 1) +
    geom_point(data = subset(plot1_data, Type == "Observed"), 
               aes(fill = Type), color = "black", shape = 21, size = 3, alpha = 1) +
    scale_color_manual(values = c("Simulated" = "green3")) +
    scale_fill_manual(values = c("Observed" = "gold")) +
    theme_minimal() +
    labs(title = paste0("Scenario concomitent_div (Iteration ", j, ")"),
         x = paste0("PC1 (", round(summary(pca1)$importance[2, 1] * 100, 1), "%)"),
         y = paste0("PC2 (", round(summary(pca1)$importance[2, 2] * 100, 1), "%)")) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  
  pca2 <- prcomp(scen2_data, scale. = TRUE, center = TRUE)
  obs_proj2 <- predict(pca2, newdata = real_new)
  
  pca2_df <- data.frame(PC1 = pca2$x[, 1], PC2 = pca2$x[, 2], Type = "Simulated")
  obs2_df <- data.frame(PC1 = obs_proj2[, 1], PC2 = obs_proj2[, 2], Type = "Observed")
  plot2_data <- rbind(pca2_df, obs2_df)
  
  p2 <- ggplot(plot2_data, aes(x = PC1, y = PC2)) +
    geom_point(data = subset(plot2_data, Type == "Simulated"), 
               aes(color = Type), alpha = 0.6, size = 1) +
    geom_point(data = subset(plot2_data, Type == "Observed"), 
               aes(fill = Type), color = "black", shape = 21, size = 3, alpha = 1) +
    scale_color_manual(values = c("Simulated" = "green3")) +
    scale_fill_manual(values = c("Observed" = "gold")) +
    theme_minimal() +
    labs(title = paste0("Scenario recent_div (Iteration ", j, ")"),
         x = paste0("PC1 (", round(summary(pca2)$importance[2, 1] * 100, 1), "%)"),
         y = paste0("PC2 (", round(summary(pca2)$importance[2, 2] * 100, 1), "%)")) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  
  pca3 <- prcomp(scen3_data, scale. = TRUE, center = TRUE)
  obs_proj3 <- predict(pca3, newdata = real_new)
  
  pca3_df <- data.frame(PC1 = pca3$x[, 1], PC2 = pca3$x[, 2], Type = "Simulated")
  obs3_df <- data.frame(PC1 = obs_proj3[, 1], PC2 = obs_proj3[, 2], Type = "Observed")
  plot3_data <- rbind(pca3_df, obs3_df)
  
  p3 <- ggplot(plot3_data, aes(x = PC1, y = PC2)) +
    geom_point(data = subset(plot3_data, Type == "Simulated"), 
               aes(color = Type), alpha = 0.6, size = 1) +
    geom_point(data = subset(plot3_data, Type == "Observed"), 
               aes(fill = Type), color = "black", shape = 21, size = 3, alpha = 1) +
    scale_color_manual(values = c("Simulated" = "green3")) +
    scale_fill_manual(values = c("Observed" = "gold")) +
    theme_minimal() +
    labs(title = paste0("Scenario old_div (Iteration ", j, ")"),
         x = paste0("PC1 (", round(summary(pca3)$importance[2, 1] * 100, 1), "%)"),
         y = paste0("PC2 (", round(summary(pca3)$importance[2, 2] * 100, 1), "%)")) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  
  pdf(paste0("PCA_plots/PCA_iteration_", j, ".pdf"), width = 15, height = 5)
  grid.arrange(p1, p2, p3, ncol = 3)
  dev.off()
  
  cat(paste0("  -> PCA plots saved: PCA_plots/PCA_iteration_", j, ".pdf\n"))
    
  # ========================================================================
  # 10 repetitions of subsamples 10,000 simulations 
  # ========================================================================
  
  confusion_matrices_repeat <- list()
  prior_errors_repeat <- numeric(10)
  model_selection_results_votes1_repeat <- list()
  model_selection_results_votes2_repeat <- list()
  model_selection_results_votes3_repeat <- list()
  model_selection_results_posterior_repeat <- list()
  
  selected_model <- list()
  selected_posterior_prob <- list()
  
  for (i in 1:10) {
    
    model_RF_i <- abcrf(formula = model ~ .,
                        data = ref_table, 
                        ntree = 1000, 
                        paral = TRUE, 
                        lda = FALSE)
    
    confusion_matrices_repeat[[i]] <- model_RF_i$model.rf$confusion.matrix
    prior_errors_repeat[i] <- model_RF_i$prior.err
    
    model_selection_result_RF <- predict(object = model_RF_i,
                                         obs = real_data,
                                         training = ref_table,
                                         ntree = 1000,
                                         paral = TRUE,
                                         paral.predict = TRUE)
    
    selected_model[[i]] <- levels(model)[which.max(model_selection_result_RF$vote)]
    selected_posterior_prob[[i]] <- model_selection_result_RF$post.prob
    
    model_selection_results_votes1_repeat[[i]] <- model_selection_result_RF$vote[1]
    model_selection_results_votes2_repeat[[i]] <- model_selection_result_RF$vote[2]
    model_selection_results_votes3_repeat[[i]] <- model_selection_result_RF$vote[3]
    model_selection_results_posterior_repeat[[i]] <- model_selection_result_RF$post.prob
  }
  
  # Calculate means and standard deviations
  global_error_rate <- mean(prior_errors_repeat)
  sd_prior_error <- sd(prior_errors_repeat)
  mean_votes1 <- mean(unlist(model_selection_results_votes1_repeat))
  mean_votes2 <- mean(unlist(model_selection_results_votes2_repeat))
  mean_votes3 <- mean(unlist(model_selection_results_votes3_repeat))
  sd_votes1 <- sd(unlist(model_selection_results_votes1_repeat))
  sd_votes2 <- sd(unlist(model_selection_results_votes2_repeat))
  sd_votes3 <- sd(unlist(model_selection_results_votes3_repeat))
  selected_model_final <- unlist(selected_model)
  selected_posterior_prob <- unlist(selected_posterior_prob)
  
  df_model <- data.frame(model = selected_model_final, post_prob = selected_posterior_prob)
  
  posterior_probabilities <- unlist(model_selection_results_posterior_repeat)
  
  max_votes <- max(mean_votes1, mean_votes2, mean_votes3)
  
  if (max_votes == mean_votes1) {
    posterior_probability <- mean(df_model$post_prob[df_model$model == "concomitent_div"])
    sd_posterior_probability <- sd(df_model$post_prob[df_model$model == "concomitent_div"])
    selected_model <- "concomitent_div"
  } else if (max_votes == mean_votes2) {
    posterior_probability <- mean(df_model$post_prob[df_model$model == "recent_div"])
    sd_posterior_probability <- sd(df_model$post_prob[df_model$model == "recent_div"])
    selected_model <- "recent_div"
  } else {
    posterior_probability <- mean(df_model$post_prob[df_model$model == "old_div"])
    sd_posterior_probability <- sd(df_model$post_prob[df_model$model == "old_div"])
    selected_model <- "old_div"
  }
  
  local_error_rate <- 1 - posterior_probability
  
  selected_model <- case_when(
    max_votes == mean_votes1 ~ "concomitent_div",
    max_votes == mean_votes2 ~ "recent_div", 
    max_votes == mean_votes3 ~ "old_div"
  )
  
  df_final <- rbind(df_final, data.frame(
    subset_data = j,
    global_error_rate = global_error_rate,
    sd_prior_error = sd_prior_error,
    mean_votes1 = mean_votes1,
    mean_votes2 = mean_votes2,
    mean_votes3 = mean_votes3,
    posterior_probability = posterior_probability,
    local_error_rate = local_error_rate,
    sd_votes1 = sd_votes1,
    sd_votes2 = sd_votes2,
    sd_votes3 = sd_votes3,
    sd_posterior_probability = sd_posterior_probability,
    selected_model = selected_model
  ))
}

print(df_final)

write.xlsx(df_final, "summary_results.xlsx", overwrite = TRUE)

cat("The Excel file 'summary_results.xlsx' has been successfully created!\n")
cat(paste0("PCA plots have been saved in the folder 'PCA_plots/'\n"))
