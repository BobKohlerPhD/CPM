library(dplyr)
library(purrr)
library(readr)
library(R.matlab)
library(jsonlite)
library(ppcor)
library(glmnet)
library(parallel)

#~~~~~~~~~~~~~~~~~~~Example Usage~~~~~~~~~~~~~~~~~~~#
# run_cpm_pipeline("config.json",
#                  "/path/to/mat_directory",
#                  "your_data.mat",
#                  "/path/to/output_directory")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~Helper Functions~~~~#
# Generate symmetric matrices from upper‐triangular vector
symm_mat <- function(vals, n) {
  m <- matrix(NA_real_, n, n)
  m[upper.tri(m)] <- vals
  m <- m + t(m)
  diag(m) <- NA_real_
  m
}

# Function to select type of model (e.g., linear, logistic, ridge)
pred_fun <- function(est, sums) {
  if (inherits(est, "cv.glmnet")) {
    as.numeric(predict(est, newx = matrix(sums, ncol = 1), s = "lambda.min"))
  } else {
    predict(est, data.frame(x = sums))
  }
}

#~~~~Determines what to calculate at each fold~~~~#
train_cpm <- function(train_mat, train_behav, n_nodes,
                      p_thresh = 0.05, mode = "linear", covariates = NULL) {
  if (!is.null(covariates)) {
    C    <- if (is.data.frame(covariates)) model.matrix(~ . - 1, covariates) else covariates
    ct   <- map(train_mat, ~ ppcor::pcor.test(.x, train_behav, C))
    r_vals <- map_dbl(ct, "estimate")
    p_vals <- map_dbl(ct, "p.value")
  } else {
    ct   <- map(train_mat, ~ cor.test(.x, train_behav))
    r_vals <- map_dbl(ct, "estimate")
    p_vals <- map_dbl(ct, "p.value")
  }
  
  # Create matrices
  r_mat <- symm_mat(r_vals, n_nodes)
  p_mat <- symm_mat(p_vals, n_nodes)
  
  # Significant edges
  pos_edges <- (r_mat > 0) & (p_mat < p_thresh)
  neg_edges <- (r_mat < 0) & (p_mat < p_thresh)
  
  # Sum edges
  ut      <- upper.tri(r_mat)
  pos_sum <- colSums(train_mat[pos_edges[ut], , drop = FALSE])
  neg_sum <- colSums(train_mat[neg_edges[ut], , drop = FALSE])
  both    <- pos_sum - neg_sum
  
  # Fit model based on the method chosen above 
  fit_one <- function(x, y) {
    if (mode == "ridge") {
      cv.glmnet(x = as.matrix(x), y = y, alpha = 0)
    } else if (mode == "linear") {
      lm(y ~ x)
    } else if (mode == "logistic") {
      glm(y ~ x, family = binomial())
    } else {
      stop("Unknown mode")
    }
  }
  
  list(    pos_est   = if (sum(pos_edges) > 0)     fit_one(pos_sum, train_behav) else NULL,
           neg_est   = if (sum(neg_edges) > 0)     fit_one(neg_sum, train_behav) else NULL,
           both_est  = if ((sum(pos_edges)|sum(neg_edges))>0) fit_one(both, train_behav) else NULL,
           pos_edges = pos_edges,
           neg_edges = neg_edges)

}

#~~~~CPM with K-FOLD CV~~~~#
kfold_cpm <- function(X, y, 
                      k = 5,
                      p_thresh = 0.05,
                      zscore = FALSE, 
                      mode = "linear", 
                      covariates = NULL) {
  n_sub   <- dim(X)[3]
  n_nodes <- dim(X)[1]
  ut      <- upper.tri(matrix(0, n_nodes, n_nodes))
  
  # flatten upper‐triangular edges
  all_edges <- apply(X, 3, function(m) m[ut])
  
  inds <- sample(n_sub)
  sz   <- floor(n_sub / k)
  out  <- list(
    y_pos   = numeric(n_sub),
    y_neg   = numeric(n_sub),
    y_both  = numeric(n_sub),
    fit_p   = vector("list", k),
    fit_n   = vector("list", k),
    fit_b   = vector("list", k),
    edges_p = vector("list", k),
    edges_n = vector("list", k)
  )
  
  for (fold in seq_len(k)) {
    test_i  <- if (fold < k) inds[((fold-1)*sz+1):(fold*sz)] else inds[((fold-1)*sz+1):n_sub]
    train_i <- setdiff(inds, test_i)
    
    tr_mat <- all_edges[, train_i, drop = FALSE]
    te_mat <- all_edges[, test_i,  drop = FALSE]
    tr_y   <- y[train_i]
    
    if (zscore) {
      tr_mat <- scale(tr_mat)
      te_mat <- scale(te_mat)
    }
    
    cpm <- train_cpm(tr_mat, tr_y, n_nodes, p_thresh, mode, covariates)
    
    ps <- rowSums(te_mat[cpm$pos_edges[ut], , drop = FALSE])
    ns <- rowSums(te_mat[cpm$neg_edges[ut], , drop = FALSE])
    bs <- ps - ns
    
    out$y_pos[test_i]  <- pred_fun(cpm$pos_est, ps)
    out$y_neg[test_i]  <- pred_fun(cpm$neg_est, ns)
    out$y_both[test_i] <- pred_fun(cpm$both_est, bs)
    
    out$fit_p[[fold]]   <- cpm$pos_est
    out$fit_n[[fold]]   <- cpm$neg_est
    out$fit_b[[fold]]   <- cpm$both_est
    out$edges_p[[fold]] <- cpm$pos_edges
    out$edges_n[[fold]] <- cpm$neg_edges
  }
  
  out
}

#~~~~Full Function~~~~#
run_cpm_pipeline <- function(k, p_thresh, zscore, mode, num_iter,
                             mat_dir, mat_name,
                             outcome_dir, outcome_name, 
                             output_dir) {
  
  matrices  <- readMat(file.path(mat_dir, mat_name))  # load matrices from .mat file 
  outcome <- read.delim(file.path(outcome_dir, outcome_name)) # load outcome of interest. file should have two columns with labels: 'subject' and 'y'
  connectivity    <- matrices$x
  outcome <- outcome$y
  subj <- outcome$subject
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  writeLines(subj, file.path(output_dir, "subjects.txt"))
  
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, c("kfold_cpm","train_cpm","X","y",
                      "symm_mat","ppcor","glmnet","pred_fun"))
  
  parLapply(cl, seq_len(num_iter), function(i) {
    res <- kfold_cpm(connectivity, outcome,
                     k        = k,
                     p_thresh = p_thresh,
                     zscore   = zscore,
                     mode     = mode)
    readr::write_csv(
      tibble(
        iter    = i,
        y_pos   = res$y_pos,
        y_neg   = res$y_neg,
        y_both  = res$y_both,
        y_true  = y
      ),
      file.path(output_dir, sprintf("pred_iter%02d.csv", i))
    )
  })
  
  stopCluster(cl)
  
  # Combine predictions from all iterations intosingle csv  
  preds <- list.files(output_dir, "pred_iter.*\\.csv$", full.names = TRUE) %>%
    map_dfr(read_csv)
  write_csv(preds, file.path(output_dir, "all_predictions.csv"))
  
  message("CPM RUN IS FINISHED AND FILES ARE SAVED")
}

#~~~~Example Launch~~~~#
run_cpm_pipeline(
  k        = 5,
  p_thresh = 0.01,
  zscore   = TRUE,
  mode     = "ridge",
  num_iter = 100,
  mat_dir  = "/path/to/mat",
  mat_name = "connect.mat",
  outcome_dir = "/path/to/outcome",
  outcome_name = "outcome.mat",
  output_dir  = "/path/to/output"
)