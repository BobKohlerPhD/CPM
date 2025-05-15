
#~~~~~~~~~~This script performs the following~~~~~~~~~~~#
# 1. Reads in CPM model output files for subsampling runs 
# 2. Compute performance metrics
# 3. Creates visualizations of performance
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(tidyverse)
library(purrr)
library(readr)
library(tidyr)
library(readxl)
library(writexl)
options(scipen = 999)

read_path <- "/Users/bobkohler/Desktop/hormone_cpm/hormone_cpm_output/cpm_output_pds/female_pds_rest/"

subsample_paths <- c('subsample1') # Example: subsample_paths <- c("/path/to/subsample1", "/path/to/subsample2", etc...)

#~~~~CPM Parameters that were used~~~~#
regression_type <- "linear"

n_repeats    <- 100
n_iterations <- 1000
lst_of_i     <- 1:(n_repeats + n_iterations)



#~~~~Loop through output folders for y_prediction csv files across subsamples~~~~#
compute_subsample_metrics <- function(pred_dir, regression_type, n_repeats, n_iterations) {
  lst_of_i <- 1:(n_repeats + n_iterations)
  
  r_pos  <- numeric()
  r_neg  <- numeric()
  r_both <- numeric()
  
  for (i in lst_of_i) {
    file_path <- file.path(pred_dir, sprintf("y_prediction_iter%d.csv", i))
    if (file.exists(file_path)) {
      df_tmp <- read.csv(file_path)
      r <- list(
        y_actual = as.numeric(df_tmp$y_actual),
        y_pos    = as.numeric(df_tmp$y_pred_pos),
        y_neg    = as.numeric(df_tmp$y_pred_neg),
        y_both   = as.numeric(df_tmp$y_pred_both))
      if (regression_type == "linear") {
        r_pos  <- c(r_pos,  cor(r$y_pos,  r$y_actual, method = "spearman", use = "complete.obs"))
        r_neg  <- c(r_neg,  cor(r$y_neg,  r$y_actual, method = "spearman", use = "complete.obs"))
        r_both <- c(r_both, cor(r$y_both, r$y_actual, method = "spearman", use = "complete.obs"))
      } else if (regression_type == "logistic") {
        accuracy_pos  <- mean(r$y_pos == r$y_actual, na.rm = TRUE)
        accuracy_neg  <- mean(r$y_neg == r$y_actual, na.rm = TRUE)
        accuracy_both <- mean(r$y_both == r$y_actual, na.rm = TRUE)
        r_pos  <- c(r_pos, accuracy_pos)
        r_neg  <- c(r_neg, accuracy_neg)
        r_both <- c(r_both, accuracy_both)
      }
    } else {
      warning(sprintf("File not found: %s", file_path))
    }
  }
  
  #~~~~Separate true (first n_repeats) and null (remaining iterations) performance~~~~#
  r_pos_true  <- r_pos[1:n_repeats]
  r_neg_true  <- r_neg[1:n_repeats]
  r_both_true <- r_both[1:n_repeats]
  r_pos_null  <- r_pos[(n_repeats + 1):length(r_pos)]
  r_neg_null  <- r_neg[(n_repeats + 1):length(r_neg)]
  r_both_null <- r_both[(n_repeats + 1):length(r_both)]
  
  #~~~~Calculate Model Performance~~~~#  
  mean_r_pos  <- mean(r_pos_true, na.rm = TRUE)
  mean_r_neg  <- mean(r_neg_true, na.rm = TRUE)
  mean_r_both <- mean(r_both_true, na.rm = TRUE)
  
  #~~~~Compute p-values~~~~#
  p_value_one_tail <- function(null, true_mean) {
    if (true_mean >= 0) {
      mean(null >= true_mean, na.rm = TRUE)
    } else {
      warning("Metric is not positive.")
      mean(null <= true_mean, na.rm = TRUE)
    }
  }
  p_value_two_tail <- function(null, true_mean) {
    mean(abs(null) >= abs(true_mean), na.rm = TRUE)
  }
  
  p_pos_onetail  <- p_value_one_tail(r_pos_null,  mean_r_pos)
  p_neg_onetail  <- p_value_one_tail(r_neg_null,  mean_r_neg)
  p_both_onetail <- p_value_one_tail(r_both_null, mean_r_both)
  p_pos_twotail  <- p_value_two_tail(r_pos_null,  mean_r_pos)
  p_neg_twotail  <- p_value_two_tail(r_neg_null,  mean_r_neg)
  p_both_twotail <- p_value_two_tail(r_both_null, mean_r_both)
  
  
  list(
    mean_r_pos = mean_r_pos,
    mean_r_neg = mean_r_neg,
    mean_r_both = mean_r_both,
    p_pos_onetail = p_pos_onetail,
    p_neg_onetail = p_neg_onetail,
    p_both_onetail = p_both_onetail,
    p_pos_twotail = p_pos_twotail,
    p_neg_twotail = p_neg_twotail,
    p_both_twotail = p_both_twotail
  )
}


#~~~~Performance Across Subsamples~~~~#
if (length(subsample_paths) > 0) {
  # Compute for each subsample
  metrics_list <- lapply(subsample_paths, function(path) {
    compute_subsample_metrics(path, regression_type, n_repeats, n_iterations)
  })
  
  # Average across subsamples
  overall_mean_r_pos  <- mean(sapply(metrics_list, function(x) x$mean_r_pos))
  overall_mean_r_neg  <- mean(sapply(metrics_list, function(x) x$mean_r_neg))
  overall_mean_r_both <- mean(sapply(metrics_list, function(x) x$mean_r_both))
  
  overall_p_pos_onetail  <- mean(sapply(metrics_list, function(x) x$p_pos_onetail))
  overall_p_neg_onetail  <- mean(sapply(metrics_list, function(x) x$p_neg_onetail))
  overall_p_both_onetail <- mean(sapply(metrics_list, function(x) x$p_both_onetail))
  
  overall_p_pos_twotail  <- mean(sapply(metrics_list, function(x) x$p_pos_twotail))
  overall_p_neg_twotail  <- mean(sapply(metrics_list, function(x) x$p_neg_twotail))
  overall_p_both_twotail <- mean(sapply(metrics_list, function(x) x$p_both_twotail))

  # Summary table
  df_summary <- tibble(
    Metric = c("Mean r Positive", "Mean r Negative", "Mean r Both"),
    Value = c(overall_mean_r_pos, overall_mean_r_neg, overall_mean_r_both),
    `p-value (One-Tail)` = c(overall_p_pos_onetail, overall_p_neg_onetail, overall_p_both_onetail),
    `p-value (Two-Tail)` = c(overall_p_pos_twotail, overall_p_neg_twotail, overall_p_both_twotail))
  print(df_summary)
  
} else {
  #~~~~Old code for obtaining performance of a run without subsampling~~~~#
  metrics <- compute_subsample_metrics(read_path, regression_type, n_repeats, n_iterations)

  df_summary <- tibble(
    `Null Model Iterations` = c((n_iterations), (n_iterations), (n_iterations)),
    `True Metric (Mean)`     = c(metrics$mean_r_pos, metrics$mean_r_neg, metrics$mean_r_both),
    `p-value (One-Tail)`     = c(metrics$p_pos_onetail, metrics$p_neg_onetail, metrics$p_both_onetail),
    `p-value (Two-Tail)`     = c(metrics$p_pos_twotail, metrics$p_neg_twotail, metrics$p_both_twotail))
  print(df_summary)
}
