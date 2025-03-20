# =============================================================================
# CPM Performance Analysis Script
# =============================================================================
# This script reads in CPM prediction output files from each iteration,
# computes performance metrics (Spearman correlation for linear regression,
# accuracy for logistic regression), and calculates p-values based on null
# permutations. A summary table of performance metrics is produced.
# =============================================================================

library(tidyverse)
library(purrr)
library(readr)
library(tidyr)
library(readxl)
library(writexl)

# ---------------------------
# Global Settings and Parameters
# ---------------------------
read_path      <- "/Users/bobkohler/Desktop/hormone_cpm/hormone_cpm_output/cpm_output_pds/female_pds_rest/"
options(scipen = 999)

# Set regression type ("linear" or "logistic")
regression_type <- "linear"

# Number of true CPM repeats and null (permutation) iterations
n_repeats    <- 100
n_iterations <- 1000
lst_of_i     <- 1:(n_repeats + n_iterations)

# Initialize vectors to store performance metrics
r_pos  <- numeric()
r_neg  <- numeric()
r_both <- numeric()

# ---------------------------
# Read Prediction Data from Iterations
# ---------------------------
data_list <- vector("list", length(lst_of_i))
for (i in lst_of_i) {
  file_path <- file.path(read_path, sprintf("y_prediction_iter%d.csv", i))
  if (file.exists(file_path)) {
    df_tmp <- read.csv(file_path)
    data_list[[i]] <- list(
      y_actual = as.numeric(df_tmp$y_actual),
      y_pos    = as.numeric(df_tmp$y_pred_pos),
      y_neg    = as.numeric(df_tmp$y_pred_neg),
      y_both   = as.numeric(df_tmp$y_pred_both)
    )
  } else {
    warning(sprintf("File not found: %s", file_path))
  }
}

# ---------------------------
# Calculate Performance Metrics
# ---------------------------
for (data in data_list) {
  if (is.null(data)) next  # Skip missing entries
  if (regression_type == "linear") {
    r_pos  <- c(r_pos,  cor(data$y_pos,  data$y_actual, method = "spearman", use = "complete.obs"))
    r_neg  <- c(r_neg,  cor(data$y_neg,  data$y_actual, method = "spearman", use = "complete.obs"))
    r_both <- c(r_both, cor(data$y_both, data$y_actual, method = "spearman", use = "complete.obs"))
  } else if (regression_type == "logistic") {
    accuracy_pos  <- mean(data$y_pos == data$y_actual, na.rm = TRUE)
    accuracy_neg  <- mean(data$y_neg == data$y_actual, na.rm = TRUE)
    accuracy_both <- mean(data$y_both == data$y_actual, na.rm = TRUE)
    r_pos  <- c(r_pos, accuracy_pos)
    r_neg  <- c(r_neg, accuracy_neg)
    r_both <- c(r_both, accuracy_both)
  }
}

# Separate metrics for true (first n_repeats) and null (remaining iterations) models
r_pos_true  <- r_pos[1:n_repeats]
r_neg_true  <- r_neg[1:n_repeats]
r_both_true <- r_both[1:n_repeats]
r_pos_null  <- r_pos[(n_repeats + 1):length(r_pos)]
r_neg_null  <- r_neg[(n_repeats + 1):length(r_neg)]
r_both_null <- r_both[(n_repeats + 1):length(r_both)]

# Print average performance for true models
metric_name <- ifelse(regression_type == "linear", "True Spearman r (Average)", "True Accuracy (Average)")
cat(sprintf("%s: Positive %.4f, Negative %.4f, Both %.4f\n", 
            metric_name,
            mean(r_pos_true, na.rm = TRUE),
            mean(r_neg_true, na.rm = TRUE),
            mean(r_both_true, na.rm = TRUE)))

# ---------------------------
# Compute p-values for Model Comparisons
# ---------------------------
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

p_pos_onetail  <- p_value_one_tail(r_pos_null,  mean(r_pos_true))
p_neg_onetail  <- p_value_one_tail(r_neg_null,  mean(r_neg_true))
p_both_onetail <- p_value_one_tail(r_both_null, mean(r_both_true))
p_pos_twotail  <- p_value_two_tail(r_pos_null,  mean(r_pos_true))
p_neg_twotail  <- p_value_two_tail(r_neg_null,  mean(r_neg_true))
p_both_twotail <- p_value_two_tail(r_both_null, mean(r_both_true))

# ---------------------------
# Create and Display Summary Table
# ---------------------------
df_p <- tibble(
  `Null Model Iterations` = c(length(r_pos_null), length(r_neg_null), length(r_both_null)
                              