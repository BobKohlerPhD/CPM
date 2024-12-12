library(tidyverse)
library(purrr)
library(readr)
library(tidyr)
library(readxl)
library(writexl)

read_path <- '/Users/bobkohler/Desktop/hormone_cpm/hormone_cpm_output/subsample_sex/subsample1_sex_logistic/'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~CHANGE REGRESSION TYPE BASED OUT CPM MODEL AND SET PARAMETERS~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
regression_type <- "logistic" #change to 'linear' or 'logistic'

n_repeats <- 100
n_iterations <- 1000
lst_of_i <- 1:(n_repeats + n_iterations)

r_pos <- c()
r_neg <- c()
r_both <- c()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~READ IN DATA FROM EACH ITERATION~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
data_list <- list()
for (i in lst_of_i) {
  file_path <- sprintf("%s/y_prediction_iter%d.csv", read_path, i)
  if (file.exists(file_path)) {
    df_tmp <- read.csv(file_path)
    # Transform and standardize columns to binary (1/0)
    data_list[[i]] <- list(
      y_actual = as.numeric(df_tmp$y_actual),
      y_pos = as.numeric(df_tmp$y_pred_pos),
      y_neg = as.numeric(df_tmp$y_pred_neg),
      y_both = as.numeric(df_tmp$y_pred_both)
    )
  } else {
    warning(sprintf("File not found: %s", file_path))
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~CALCULATE CORR/ACCURACY ACROSS ALL ITERATIONS~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
for (i in seq_along(data_list)) {
  data <- data_list[[i]]
  
  if (regression_type == "linear") {
    # Compute Spearman correlations
    r_pos <- c(r_pos, cor(data$y_pos, data$y_actual, method = "spearman", use = "complete.obs"))
    r_neg <- c(r_neg, cor(data$y_neg, data$y_actual, method = "spearman", use = "complete.obs"))
    r_both <- c(r_both, cor(data$y_both, data$y_actual, method = "spearman", use = "complete.obs"))
  } else if (regression_type == "logistic") {
    # Compute Logistic Accuracy (row-wise comparison)
    accuracy_pos <- sum(data$y_pos == data$y_actual, 1, 0)
      accuracy_neg <- sum(data$y_neg == data$y_actual, 1, 0)
        accuracy_both <- sum(data$y_both == data$y_actual, 1, 0)
          
          # Add mean row-wise accuracy to results
          r_pos <- c(r_pos, mean(accuracy_pos, na.rm = TRUE))
          r_neg <- c(r_neg, mean(accuracy_neg, na.rm = TRUE))
          r_both <- c(r_both, mean(accuracy_both, na.rm = TRUE))
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~SPLIT/PRINT TRUE AND NULL MODEL OUTPUT~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

r_pos_null <- r_pos[(n_repeats + 1):length(r_pos)]/length(data$y_actual)
r_pos_true <- r_pos[1:n_repeats]/length(data$y_actual)

r_neg_null <- r_neg[(n_repeats + 1):length(r_neg)]/length(data$y_actual)
r_neg_true <- r_neg[1:n_repeats]/length(data$y_actual)

r_both_null <- r_both[(n_repeats + 1):length(r_both)]/length(data$y_actual)
r_both_true <- r_both[1:n_repeats]/length(data$y_actual)

#~~PRINT RESULTS~~#
if (regression_type == "linear") {
  metric_name <- "True Spearman r (Average)"
} else {
  metric_name <- "True Accuracy (Average)"
}

cat(sprintf("%s: Positive %.4f, Negative %.4f, Both %.4f\n", 
            metric_name,
            mean(r_pos_true, na.rm = TRUE),
            mean(r_neg_true, na.rm = TRUE),
            mean(r_both_true, na.rm = TRUE)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~CALCULATE P VALUES BASED ON PERMUTATIONS~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
p_value_one_tail <- function(null, true_mean) {
  if (true_mean >= 0) {
    return(mean(null >= true_mean, na.rm = TRUE))
  } else {
    warning("Metric is not positive.")
    return(mean(null <= true_mean, na.rm = TRUE))
  }
}

p_pos_onetail <- p_value_one_tail(r_pos_null, mean(r_pos_true))
p_neg_onetail <- p_value_one_tail(r_neg_null, mean(r_neg_true))
p_both_onetail <- p_value_one_tail(r_both_null, mean(r_both_true))

p_value_two_tail <- function(null, true_mean) {
  abs_null <- abs(null)
  return(mean(abs_null >= abs(true_mean), na.rm = TRUE))
}

p_pos_twotail <- p_value_two_tail(r_pos_null, mean(r_pos_true))
p_neg_twotail <- p_value_two_tail(r_neg_null, mean(r_neg_true))
p_both_twotail <- p_value_two_tail(r_both_null, mean(r_both_true))

#~~~~~~~~~~~~~~~~~~~~~~~~#
#~~SUMMARY TABLE OUTPUT~~#
#~~~~~~~~~~~~~~~~~~~~~~~~#
(df_p <- tibble(
  `Null Model Iterations` = c(length(r_pos_null), length(r_neg_null), length(r_both_null)),
  `True Metric (Mean)` = c(mean(r_pos_true, na.rm = TRUE), mean(r_neg_true, na.rm = TRUE), mean(r_both_true, na.rm = TRUE)),
  `True Metric (SD)` = c(sd(r_pos_true, na.rm = TRUE), sd(r_neg_true, na.rm = TRUE), sd(r_both_true, na.rm = TRUE)),
  `p-value (One-Tail)` = c(p_pos_onetail, p_neg_onetail, p_both_onetail),
  `p-value (Two-Tail)` = c(p_pos_twotail, p_neg_twotail, p_both_twotail)))

(df_boxplot <- tibble(
  `Null Model Size` = c(length(r_pos_null), length(r_neg_null), length(r_both_null)),
  `True Model Size` = c(length(r_pos_true), length(r_neg_true), length(r_both_true)),
  `Metric (Mean)` = c(mean(r_pos_true, na.rm = TRUE), mean(r_neg_true, na.rm = TRUE), mean(r_both_true, na.rm = TRUE)),
  `p-value (One-Tail)` = c(p_pos_onetail, p_neg_onetail, p_both_onetail),
  `Metric (Median)` = c(median(r_pos_true, na.rm = TRUE), median(r_neg_true, na.rm = TRUE), median(r_both_true, na.rm = TRUE))))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~CREATE BOXPLOTS FOR ACCURACY OVER ITERATIONS~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
level_order <- c('Positive', 'Negative', 'Both') 
create_plot <- function(data, x_var, y_var, edge_var, level_order, fill_values, color_values) {
  ggplot(data, aes(x = factor(!!rlang::sym(x_var), levels = level_order), 
                   y = !!rlang::sym(y_var), color = !!rlang::sym(edge_var), fill = !!rlang::sym(edge_var))) +
    geom_boxplot(aes(fill = !!rlang::sym(edge_var),
                     fill = after_scale(colorspace::lighten(fill, .2))),
                 color = "black",
                 width = .55,
                 fatten = .75,
                 linewidth = 1.25,
                 outlier.shape = NA)+
    geom_point(position = position_jitter(width = .2, seed = 0), size = 4, alpha = .4) +
    geom_point(position = position_jitter(width = .2, seed = 0), size = 4, stroke = .4, shape = 1) +
    scale_fill_manual(values = fill_values) +
    scale_color_manual(values = color_values) +
    theme_minimal() +
    labs(x = "", y = "Accuracy" , title = "Average Cross-Validated Performance by Edge Type")+
    theme(text = element_text(family = "Arial"),
          plot.title = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 22, face = "bold"),
          axis.line = element_line(linewidth = 1.25),
          axis.ticks.length = unit(.25, "cm"),
          axis.text.y = element_text(size = 18, color = "black"),
          axis.text.x = element_text(size = 22, color = "black", face = "bold"),
          legend.text = element_text(size = 0, color = "black", face = "bold"),
          legend.title = element_text(size = 0, color = "black", face = "bold"),
          legend.position = "none")
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~PLOT TRUE CORRELATIONS FOR REPEATS~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
r_true_df <- data.frame(Positive = r_pos_true, Negative = r_neg_true, Both = r_both_true) %>% 
  pivot_longer(
    cols = everything(),
    names_to = "edge",
    values_to = "r_value",
    values_drop_na = FALSE)

(r_true_boxplot <- 
  create_plot(r_true_df, 
              x_var = "edge", 
              y_var = "r_value", 
              edge_var = "edge", 
              level_order = level_order, 
              fill_values = c('Positive' = 'red', 'Negative' = 'blue', 'Both' = 'purple'), 
              color_values = c('Positive' = 'red', 'Negative' = 'blue', 'Both' = 'purple')))


