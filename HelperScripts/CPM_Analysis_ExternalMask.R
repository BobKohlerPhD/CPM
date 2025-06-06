# This script applies an existing connectivity mask
# (pos_edges, neg_edges) to a new set of matrices 
# and performs correlations appropriate analysis based on datatype 



# x_array: [n_nodes × n_nodes × n_subj] connectivity matrix 
# y_vector: outcome of interest for new sample that mask will be applied to
# pos_edges: [n_nodes × n_nodes] positive mask from original CPM 
# neg_edges: [n_nodes × n_nodes] negative mask from original CPM 
# outcome_type: "continuous", "binary" 

library(dplyr)
library(tibble)
library(pROC)


# Load in new matrices - not tested on .mat files 
new_matrices <- readRDS("path/to/new/matrices/for/applying/external/mask")


# Data Paths 
filtered_outcome_path <- 'path/to/outcome/filtered/for/good/matrices'
pos_raw_path <- 'path/to/pos/mask.csv'
neg_raw_path <- 'path/to/neg/mask.csv'


# Path and name of file for results 
output_results_path <- 'path/to/where/you/want/results'
output_results_filename <- "results_externalmask.rds"


subject_ids <- names(new_matrices)
n_subj  <- length(new_matrices)
n_nodes <- nrow(new_matrices[[1]])
x_array <- array(NA_real_, dim = c(n_nodes, n_nodes, n_subj))

for (i in seq_len(n_subj)) {
  mat_i <- new_matrices[[i]]
  stopifnot(is.matrix(mat_i),
            nrow(mat_i) == n_nodes,
            ncol(mat_i) == n_nodes)
  x_array[ , , i] <- mat_i
}


outcomes_df <- read.csv(filtered_outcome_path) %>%
  mutate(subject_id = src_subject_id) %>%
  mutate(outcome = any.alcohol) %>%
  select(-X, -src_subject_id)

stopifnot(all(c("subject_id", "outcome") %in% names(outcomes_df)))

outcome_vec <- outcomes_df$outcome
names(outcome_vec) <- outcomes_df$subject_id

# Check every subject has corresponding outcome
stopifnot(all(subject_ids %in% names(outcome_vec)))

y_vec <- outcome_vec[subject_ids]
stopifnot(length(y_vec) == n_subj)


pos_raw = read.table(pos_raw_path,
                       header = FALSE,
                       stringsAsFactors = FALSE)

neg_raw = read.table(neg_raw_path, 
                       header = FALSE, 
                       stringsAsFactors = FALSE) 

pos_mat <- as.matrix(pos_raw) 
neg_mat <- as.matrix(neg_raw)

# Convert masks to True/False instead of 0/1 
pos_edges <- pos_mat != 0
neg_edges <- neg_mat != 0

n_nodes <- 268 
stopifnot(is.matrix(pos_edges), all(dim(pos_edges) == c(n_nodes, n_nodes)))
stopifnot(is.matrix(neg_edges), all(dim(neg_edges) == c(n_nodes, n_nodes)))



#~~~~Function for applying external mask and calculating relevant statistics~~~~#
apply_external_mask <- function(x_array, y_vec, pos_edges, neg_edges,
                                       outcome_type = c("continuous", "binary")) {
  outcome_type <- match.arg(outcome_type)
  
  # Checks dimensions of matrix files to be safe 
  stopifnot(is.array(x_array), length(dim(x_array)) == 3)
  n_nodes <- dim(x_array)[1]
  n_subj  <- dim(x_array)[3]
  
  
  stopifnot(length(y_vec) == n_subj)
  stopifnot(is.matrix(pos_edges), all(dim(pos_edges) == c(n_nodes, n_nodes)))
  stopifnot(is.matrix(neg_edges), all(dim(neg_edges) == c(n_nodes, n_nodes)))
  
  # Upper-triangular edges 
  upper_tri  <- upper.tri(matrix(0, n_nodes, n_nodes))
  edge_mat   <- apply(x_array, 3, function(m) m[upper_tri])  
  mask_p_vec <- pos_edges[upper_tri]
  mask_n_vec <- neg_edges[upper_tri]
  
  # Calculate sum strength
  sum_pos  <- colSums(edge_mat[mask_p_vec, , drop = FALSE])
  sum_neg  <- colSums(edge_mat[mask_n_vec, , drop = FALSE])
  sum_both <- sum_pos - sum_neg
  
  # DF with sum strength and outcome for final analyses 
  df_strength <- tibble(
    subject = seq_len(n_subj),
    pos_strength  = sum_pos,
    neg_strength  = sum_neg,
    both_strength = sum_both,
    outcome       = y_vec)

  
  # Continuous outcome uses pearson 
  if (outcome_type == "continuous") { 
    cor_pos  <- cor.test(sum_pos,  y_vec, method = "pearson") #can change to spearman depending on outcome distribution 
    cor_neg  <- cor.test(sum_neg,  y_vec, method = "pearson")
    cor_both <- cor.test(sum_both, y_vec, method = "pearson")
    
    return(list(
      strength = df_strength,
      result_pos  = cor_pos,
      result_neg  = cor_neg,
      result_both = cor_both
    ))
  }
  
  # Binary outcome uses logistic regression & t-test  
  df_strength <- df_strength %>%
    mutate(outcome = factor(outcome)) 
  
  # Logistic 
  glm_pos  <- glm(outcome ~ pos_strength,  data = df_strength, family = binomial(link = "logit"))
  glm_neg  <- glm(outcome ~ neg_strength,  data = df_strength, family = binomial(link = "logit"))
  glm_both <- glm(outcome ~ both_strength, data = df_strength, family = binomial(link = "logit"))
  
  # p-values 
  p_pos  <- summary(glm_pos)$coefficients["pos_strength", "Pr(>|z|)"]
  p_neg  <- summary(glm_neg)$coefficients["neg_strength", "Pr(>|z|)"]
  p_both <- summary(glm_both)$coefficients["both_strength", "Pr(>|z|)"]
  
  # AUC 
  roc_pos  <- roc(df_strength$outcome, predict(glm_pos,  type = "response"))
  roc_neg  <- roc(df_strength$outcome, predict(glm_neg,  type = "response"))
  roc_both <- roc(df_strength$outcome, predict(glm_both, type = "response"))
  
  # t-test
  t_pos  <- t.test(pos_strength  ~ outcome, data = df_strength)
  t_neg  <- t.test(neg_strength  ~ outcome, data = df_strength)
  t_both <- t.test(both_strength ~ outcome, data = df_strength)
  
  return(list(
    strength    = df_strength,
    glm_models  = list(pos = glm_pos, neg = glm_neg, both = glm_both),
    p_values    = tibble(pos = p_pos, neg = p_neg, both = p_both),
    auc         = tibble(pos = auc(roc_pos), neg = auc(roc_neg), both = auc(roc_both)),
    t_tests     = list(pos = t_pos, neg = t_neg, both = t_both)
  ))
}


results_externalmask <- apply_external_mask(x_array      = x_array,
                                            y_vec        = y_vec,
                                            pos_edges    = pos_edges,
                                            neg_edges    = neg_edges,
                                            outcome_type = "binary")



saveRDS(results_externalmask,
        file.path(output_results_path, output_results_filename))


