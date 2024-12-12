library(caret)
library(dplyr)
library(purrr)

#~~~~~~~~~~~~~~~~#
#~~CPM FUNCTION~~#
#~~~~~~~~~~~~~~~~#

train_cpm <- function(train_mat, train_behav, num_nodes, p_thresh = 0.05, mode = "linear") {
  #~~COMPUTE CORRELATIONS~~#
  corr_train <- map(train_mat, ~ cor.test(.x, train_behav, method = "pearson"))
  r_lst <- sapply(corr_train, `[[`, "estimate")
  p_lst <- sapply(corr_train, `[[`, "p.value")
  
  #~~MATRIX SETUP~~#
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
  
  #~~CHECK SYMMETRY OF MATRICES~~#
  if (!isSymmetric(r_mat) || !isSymmetric(p_mat)) {
    stop("ERROR: r_mat or p_mat is not symmetric. Please check your data.")
  }
  
  #~~IDENTIFY SIGNIFICANT EDGES~~#
  pos_edges <- (r_mat > 0) & (p_mat < p_thresh)
  neg_edges <- (r_mat < 0) & (p_mat < p_thresh)
  
  #~~SUM EDGES~~#
  pos_sum <- colSums(train_mat[pos_edges[upper_tri], , drop = FALSE])
  neg_sum <- colSums(train_mat[neg_edges[upper_tri], , drop = FALSE])
  both <- pos_sum - neg_sum
  
  #~~MODEL TRAINING~~#
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
  both_estimator <- if (sum(pos_edges | neg_edges) > 0) train_model(train_behav, both, mode) else NA
  
  list(pos_estimator = pos_estimator, 
       neg_estimator = neg_estimator, 
       both_estimator = both_estimator, 
       pos_edges = pos_edges, 
       neg_edges = neg_edges)
}

#~~KFOLD CV FUNCTION~~#
kfold_cpm <- function(x, y, k, p_thresh = 0.05, zscore = FALSE, mode = "linear") {
  num_subs <- dim(x)[3]
  num_nodes <- dim(x)[1]
  upper_tri <- upper.tri(matrix(0, num_nodes, num_nodes))
  
  #~~MATRIX FLATTEN~~#
  all_edges <- apply(x, 3, function(mat) mat[upper_tri])
  
  rand_inds <- sample(seq_len(num_subs))
  sample_size <- floor(num_subs / k)
  
  results <- list(
    y_pred_pos = numeric(num_subs),
    y_pred_neg = numeric(num_subs),
    y_pred_both = numeric(num_subs),
    fit_p = vector("list", k),
    fit_n = vector("list", k),
    fit_b = vector("list", k),
    edges_p = vector("list", k),
    edges_n = vector("list", k)
  )
  
  for (fold in seq_len(k)) {
    test_inds <- if (fold != k) {
      rand_inds[((fold - 1) * sample_size + 1):(fold * sample_size)]
    } else {
      rand_inds[((fold - 1) * sample_size + 1):num_subs]
    }
    train_inds <- setdiff(rand_inds, test_inds)
    
    train_mats <- all_edges[, train_inds]
    train_behav <- y[train_inds]
    test_mats <- all_edges[, test_inds]
    
    if (zscore) {
      train_mats <- scale(train_mats)
      test_mats <- scale(test_mats)
    }
    
    cpm <- train_cpm(train_mats, train_behav, num_nodes, p_thresh, mode)
    pos_sum <- rowSums(test_mats[cpm$pos_edges[upper_tri], , drop = FALSE])
    neg_sum <- rowSums(test_mats[cpm$neg_edges[upper_tri], , drop = FALSE])
    both <- pos_sum - neg_sum
    
    results$y_pred_pos[test_inds] <- predict(cpm$pos_estimator, pos_sum)
    results$y_pred_neg[test_inds] <- predict(cpm$neg_estimator, neg_sum)
    results$y_pred_both[test_inds] <- predict(cpm$both_estimator, both)
    results$fit_p[[fold]] <- cpm$pos_estimator
    results$fit_n[[fold]] <- cpm$neg_estimator
    results$fit_b[[fold]] <- cpm$both_estimator
    results$edges_p[[fold]] <- cpm$pos_edges
    results$edges_n[[fold]] <- cpm$neg_edges
  }
  
  results
}

#~~RUN PARALLELIZED CPM~~#
run_cpm_thread <- function(y, iter, x, k, out_path, p_thresh = 0.05, zscore = FALSE, mode = "linear") {
  message(sprintf("Running iteration #%d", iter))
  results <- kfold_cpm(x, y, k, p_thresh, zscore, mode)
  
  saveRDS(results, file.path(out_path, sprintf("results_iter_%d.rds", iter)))
}
