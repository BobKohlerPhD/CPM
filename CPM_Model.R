library(caret)
library(dplyr)
library(purrr)
library(readr)
library(R.matlab)
library(jsonlite)
library(parallel)
library(fs)
library(stringr)
library(lubridate)
library(ppcor) 

#============================================================#
# UTILITY FUNCTIONS (CPM functions and helpers)              #
#============================================================#

#~~Model Functions~~#
train_cpm <- function(train_mat, train_behav, num_nodes, p_thresh = 0.05, mode = "linear", covariates = NULL) {
  # Identify if covariates are included
  if (!is.null(covariates)) {
    if (is.data.frame(covariates)) {
      # Dummy code categorical 
      cov_mat <- model.matrix(~ . - 1, data = covariates)
    } else {
      cov_mat <- covariates
    }
  } else {
    cov_mat <- NULL
  }
  
  # Compute correlations (if covariates are provided, use partial correlations)
  if (!is.null(cov_mat)) {
    # partial correlation  and p-values for each edge
    corr_train <- map(train_mat, ~ ppcor::pcor.test(.x, train_behav, cov_mat))
    r_lst <- sapply(corr_train, function(res) res$estimate)
    p_lst <- sapply(corr_train, function(res) res$p.value)
  } else {
    corr_train <- map(train_mat, ~ cor.test(.x, train_behav, method = "pearson"))
    r_lst <- sapply(corr_train, `[[`, "estimate")
    p_lst <- sapply(corr_train, `[[`, "p.value")
  }
  
  # Matrix Setup
  r_mat <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
  p_mat <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
  upper_tri <- upper.tri(r_mat)
  lower_tri <- lower.tri(r_mat)
  diag(r_mat) <- NA
  diag(p_mat) <- NA
  
  r_mat[upper_tri] <- r_lst
  p_mat[upper_tri] <- p_lst
  r_mat[lower_tri] <- t(r_mat)[lower_tri]
  p_mat[lower_tri] <- t(p_mat)[lower_tri]
  
  # Check Matrix Symmetry
  if (!isSymmetric(r_mat) || !isSymmetric(p_mat)) {
    stop("ERROR: r_mat or p_mat is not symmetric. Please check your data.")
  }
  
  # Identify Significant Edges
  pos_edges <- (r_mat > 0) & (p_mat < p_thresh)
  neg_edges <- (r_mat < 0) & (p_mat < p_thresh)
  
  # Sum edges (using the significant edges in the upper triangle)
  pos_sum <- colSums(train_mat[pos_edges[upper_tri], , drop = FALSE])
  neg_sum <- colSums(train_mat[neg_edges[upper_tri], , drop = FALSE])
  both <- pos_sum - neg_sum
  
  #  Train Function
  train_model <- function(y, x, mode) {
    if (mode == "ridge") {
      train(x, y, method = "ridge", tuneLength = 10)
    } else if (mode == "linear") {
      lm(y ~ x)
    } else if (mode == "logistic") {
      glm(y ~ x, family = binomial())
    } else {
      stop("Mode not implemented!")
    }
  }
  
  pos_estimator <- if (sum(pos_edges) > 0) train_model(train_behav, pos_sum, mode) else NA
  neg_estimator <- if (sum(neg_edges) > 0) train_model(train_behav, neg_sum, mode) else NA
  both_estimator <- if (sum(pos_edges | neg_edges) > 0) train_model(train
                                                                    