# Generate the figure 6

library(poolfstat)
library(tidyr)
library(dplyr)
library(ggpubr)
library(ape)
library(igraph)
library(ggplot2)
library(zoo)
library(svglite)
library(gridExtra)
library(scales)
library(grid)
library(gtable)
library(openxlsx)
options(scipen = 999)



#-----------------------------------------------------------
# Create plots for genomic regions
#-----------------------------------------------------------

# Load Baypass result
load("path/to/baypass_results.RData")
source("path/to/baypass_utils.R")
q1 <- quantile(local_score_fb_c2_SP_WP$res.local.scores$lindley, 0.99)

# Define genomic windows
position_pairs <- list(
  c(12000000, 14600000),
  c(16600000, 18800000)
)

for (z in 1:length(position_pairs)) {
  pair <- position_pairs[[z]]
  start_win <- pair[1]
  end_win <- pair[2]
  
  # Load gene annotations
  genes_Z_10 <- read.xlsx(paste0("path/to/gene_annotations/", start_win, "_", 
                                 end_win, "_genes.xlsx"))
  
  # Assign gene numbers and positions
  if (start_win == 12000000) {
    genes_Z_10$numero <- 1:nrow(genes_Z_10)
    first_num <- 1
    last_num <- 60
  } else {
    genes_Z_10$numero <- 61:(60 + nrow(genes_Z_10))
    first_num <- 61
    last_num <- 103
  }
  
  genes_Z_10$vertical_position <- rep(1:4, length.out = nrow(genes_Z_10))
  genes_Z_10$show_number <- genes_Z_10$circadian %in% c("yes", "maybe") | 
                           genes_Z_10$numero == first_num | 
                           genes_Z_10$numero == last_num
  
  # Plot 1: Gene track
  graph1 <- ggplot() +
    geom_rect(data = genes_Z_10,
              aes(xmin = start, xmax = end,
                  ymin = 1,
                  ymax = ifelse(vertical_position == 1 | vertical_position == 3, 1.07, 0.93),
                  fill = circadian), color = "black", size = 0.15) +
    geom_text(data = subset(genes_Z_10, show_number == TRUE),
              aes(x = (start + end)/2,
                  y = ifelse(vertical_position == 1, 1.15, 
                            ifelse(vertical_position == 2, 0.86,
                                  ifelse(vertical_position == 3, 1.15, 0.86))),
                  label = numero, color = circadian),
              size = 12) +
    annotate("segment", x = min(genes_Z_10$start), xend = max(genes_Z_10$end),
             y = 1, yend = 1, color = "black", size = 1) +
    scale_fill_manual(values = c("yes" = "red", "no" = "grey80", "maybe" = "blue")) +
    scale_color_manual(values = c("yes" = "red", "no" = "grey31", "maybe" = "blue")) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 34),
      axis.line = element_line(color = "white", size = 0.75),
      axis.ticks = element_line(color = "white", size = 0.75),
      axis.ticks.length = unit(0.4, "cm"),
      axis.text.y = element_text(size = 32, color = "white"),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 36, color = "white", margin = margin(r = 15)),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.position = "none"
    ) +
    labs(y = expression(F[ST] ~ "values")) +
    scale_y_continuous(limits = c(0.80, 1.20), 
                      labels = function(x) format(x, nsmall = 2)) +
    coord_cartesian(xlim = c(start_win, end_win))
  
  # Plot 2: Baypass C2 local scores
  graph2 <- ggplot(data = local_score_fb_c2_SP_WP$res.local.scores, 
                  aes(x = pos, y = lindley)) +
    geom_point(color = "grey", size = 2.5, shape = 16) +
    geom_hline(yintercept = q1, linetype = "dashed", 
              color = alpha("orange", 0.8), size = 1) +
    scale_y_continuous(limits = c(0, 60), 
                      labels = function(x) format(x, nsmall = 2)) +
    xlim(c(start_win, end_win)) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 32),
      axis.title = element_text(size = 34),
      axis.title.y = element_text(size = 34, margin = margin(r = 15)),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 0.75),
      axis.ticks = element_line(color = "black", size = 0.75),
      axis.ticks.length = unit(0.4, "cm"),
      panel.border = element_blank()
    ) +
    labs(title = "", x = "", y = expression(C[2]~"Local-score"))
  
  # Plot 3: FST scan
  load("/path/to/fst/result/fst_sliding.RData")
  quantile_95_fst <- quantile(fst_df$fst, 0.95)
  
  graph3 <- ggplot(fst_df, aes(x = pos, y = fst)) +
    geom_point(color = "grey", size = 2, alpha = 0.5) +
    geom_line(data = rolling_df, aes(x = pos, y = mean), 
             color = alpha("black", 0.6), size = 2) +
    geom_point(data = df_cat_fst_purple, aes(x = pos, y = fst, size = 7),
              color = df_cat_fst_purple$color, shape = 8) +
    geom_point(data = df_cat_fst_green, aes(x = pos, y = fst, size = 7),
              color = df_cat_fst_green$color, shape = 8) +
    geom_point(data = df_cat_fst_red, aes(x = pos, y = fst, size = 7),
              color = df_cat_fst_red$color, shape = 8) +
    scale_y_continuous(limits = c(-0.035, 1), 
                      labels = function(x) format(x, nsmall = 3)) +
    xlim(c(start_win, end_win)) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 32),
      axis.title = element_text(size = 34),
      axis.title.y = element_text(size = 34, margin = margin(r = 15)),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 0.75),
      axis.ticks = element_line(color = "black", size = 0.75),
      axis.ticks.length = unit(0.4, "cm"),
      panel.border = element_blank(),
      legend.position = "none"
    ) +
    geom_hline(yintercept = 1, color = "black") +
    labs(title = "", x = "", y = expression(F[ST] ~ "values"))
  
  # Plot 4: Nucleotide diversity (Pi)
  load("path/to/pi_data.RData")
  pi_SP$mid <- (pi_SP$window_pos_1 + pi_SP$window_pos_2) / 2
  pi_WP$mid <- (pi_WP$window_pos_1 + pi_WP$window_pos_2) / 2
  
  rolling_mean_pi <- rollmean(zoo(pi_SP$avg_pi, order.by = pi_SP$mid), 
                              k = 20, na.rm = TRUE)
  rolling_df_pi <- data.frame(pos = index(rolling_mean_pi),
                              mean = coredata(rolling_mean_pi))
  
  rolling_mean_pi_WP <- rollmean(zoo(pi_WP$avg_pi, order.by = pi_WP$mid), 
                                k = 20, na.rm = TRUE)
  rolling_df_pi_WP <- data.frame(pos = index(rolling_mean_pi_WP),
                                mean = coredata(rolling_mean_pi_WP))
  
  graph4 <- ggplot(pi_SP, aes(mid, avg_pi)) +
    xlim(c(start_win, end_win)) +
    ylim(c(0, 0.008)) +
    geom_point(color = alpha("brown1", 0.1), size = 2) +
    geom_point(data = pi_WP, aes(mid, avg_pi), 
              color = alpha("cadetblue3", 0.1), size = 2) +
    geom_line(data = rolling_df_pi, aes(x = pos, y = mean), 
             color = alpha("brown1", 0.9), size = 1.5) +
    geom_line(data = rolling_df_pi_WP, aes(x = pos, y = mean), 
             color = alpha("cadetblue3", 0.9), size = 1.5) +
    scale_color_identity() +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 32),
      axis.title = element_text(size = 34),
      axis.title.y = element_text(size = 34, margin = margin(r = 15)),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 0.75),
      axis.ticks = element_line(color = "black", size = 0.75),
      axis.ticks.length = unit(0.4, "cm"),
      panel.border = element_blank()
    ) +
    labs(title = "", x = "", y = expression(pi ~ "values"))
  
  # Plot 5: Tajima's D
  D_taj <- read.table(paste0("path/to/tajima_d/chrZ_LSP_MAF005_2000.Tajima.D"), h = T)
  D_taj_WP <- read.table(paste0("path/to/tajima_d/chrZ_LWP_MAF005_2000.Tajima.D"), h = T)
  
  rolling_mean_taj <- rollmean(zoo(D_taj$TajimaD, order.by = D_taj$BIN_START),
                               k = 75, na.rm = TRUE)
  rolling_df_taj <- data.frame(pos = index(rolling_mean_taj),
                               mean = coredata(rolling_mean_taj))
  
  rolling_mean_taj_WP <- rollmean(zoo(D_taj_WP$TajimaD, order.by = D_taj_WP$BIN_START),
                                 k = 75, na.rm = TRUE)
  rolling_df_taj_WP <- data.frame(pos = index(rolling_mean_taj_WP),
                                 mean = coredata(rolling_mean_taj_WP))
  
  graph5 <- ggplot(D_taj, aes(BIN_START, TajimaD)) +
    xlim(c(start_win, end_win)) +
    geom_point(color = alpha("brown1", 0.1), size = 2) +
    geom_point(data = D_taj_WP, aes(BIN_START, TajimaD), 
              color = alpha("cadetblue3", 0.1), size = 2) +
    geom_line(data = rolling_df_taj, aes(x = pos, y = mean), 
             color = alpha("brown1", 0.9), size = 1.5) +
    geom_line(data = rolling_df_taj_WP, aes(x = pos, y = mean), 
             color = alpha("cadetblue3", 0.9), size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", 
              color = alpha("black", 0.4), size = 2) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 3)) +
    scale_color_identity() +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 32),
      axis.title = element_text(size = 34),
      axis.title.y = element_text(size = 34, margin = margin(r = 15)),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 0.75),
      axis.ticks = element_line(color = "black", size = 0.75),
      axis.ticks.length = unit(0.4, "cm"),
      panel.border = element_blank()
    ) +
    labs(title = "", x = "", y = "Tajima's D")
  
  # Plot 6: Read depth ratio
  load("path/to/coverage_data.RData")
  
  quantile_5_couv <- quantile(result_df_SP$ratio, 0.01, na.rm = T)
  quantile_95_couv <- quantile(result_df_SP$ratio, 0.99, na.rm = T)
  
  colors_couv <- ifelse(result_df_SP$ratio <= quantile_5_couv, 
                       scales::alpha("orange", 0.6),
                       ifelse(result_df_SP$ratio >= quantile_95_couv, 
                              scales::alpha("orange", 0.6),
                              scales::alpha("grey", 0.6)))
  
  graph6 <- ggplot(result_df_SP, aes(x = V2, y = ratio)) +
    xlim(c(start_win, end_win)) +
    geom_point(aes(color = colors_couv), size = 2.5) +
    scale_y_continuous(limits = c(-3.1, 1), 
                      labels = function(x) format(x, nsmall = 3)) +
    scale_color_identity() +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 32),
      axis.title.x = element_text(size = 36, margin = margin(t = 20)),
      axis.title.y = element_text(size = 34, margin = margin(r = 15)),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 0.75),
      axis.ticks = element_line(color = "black", size = 0.75),
      axis.ticks.length = unit(0.4, "cm"),
      panel.border = element_blank()
    ) +
    labs(title = "", x = "Position on chromosome Z (in bp)", 
         y = "Ratio read depth\n difference LSP vs LWP")
  
  # Store plots for this region
  plots_for_region <- list(graph1, graph2, graph3, graph4, graph5, graph6)
  
  # Convert to gtable
  gt_list <- lapply(plots_for_region, ggplotGrob)
  
  # Find maximum y-axis width
  max_width <- do.call(grid::unit.pmax, lapply(gt_list, function(x) {
    x$widths[2:5]
  }))
  
  # Adjust all plots to this maximum width
  for (j in seq_along(gt_list)) {
    gt_list[[j]]$widths[2:5] <- max_width
  }
  
  # Assign plots to specific column
  if (z == 1) {
    plots_column1 <- gt_list
  } else {
    plots_column2 <- gt_list
  }
}

#-----------------------------------------------------------
# Create composite figure
#-----------------------------------------------------------

# Create letter plots
create_letter_plot <- function(letter) {
  ggplot() + 
    annotate("text", x = 0.5, y = 0.95, label = paste0(letter, ")"), 
             size = 18, hjust = 0.5, vjust = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1)) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

letters <- c("A", "C", "E", "G", "I", "K")
letters_B <- c("B", "D", "F", "H", "J", "L")

letter_plots_col1 <- lapply(letters[1:6], create_letter_plot)
letter_plots_col2 <- lapply(letters_B[1:6], create_letter_plot)

# Combine all plots
all_plots <- c(letter_plots_col1, plots_column1, letter_plots_col2, plots_column2)

# Layout matrix
layout_matrix <- matrix(1:24, nrow = 6, ncol = 4)

# Save composite figure
pdf(file = "Genomic_Regions_Composite.pdf",
    width = 2600/72, height = 2200/72)

grid.arrange(
  grobs = all_plots, 
  ncol = 4, 
  widths = c(0.03, 0.44, 0.03, 0.44),
  heights = c(2.25, 7, 7, 7, 7, 7),
  layout_matrix = layout_matrix
)

dev.off()
