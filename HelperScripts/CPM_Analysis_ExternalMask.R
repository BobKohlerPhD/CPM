library(dplyr)
library(tibble)
library(pROC)


# This script applies an existing connectivity mask
# (pos_edges, neg_edges) to a new set of matrices 
# and performs correlations based on datatype. AUC is also 
# generated and other stats can be easily implemented as needed 


# Data Paths for new matrices for applying external mask and outcome variable 
new_matrices_path <- '/Users/bobkohler/Desktop/ARS/ABCD_Projects/abcd_matrices_general/timepoint2/tp2_sst_filteredY3NEGEXP_matrices.rds'
filtered_outcome_path <- '/Users/bobkohler/Desktop/ARS/ABCD_Projects/alcohol_cpm/data/tp2_sst_filtered_Y3NEGEXP_outcome.csv'

# Alslows function to change these strings to a common name
outcome_var <- "aeq_positive_expectancies_ss.3_year" # name of outcome variable 
outcome_type <- "continuous" # choose continuous or binary for outcome 
subject_id_var <- "src_subject_id" # name of subject id variable 

# Path to external masks 
pos_raw_path <- '/Users/bobkohler/Desktop/ARS/yip_jama_psych/use these/pos_mat_both.csv'
neg_raw_path <- '/Users/bobkohler/Desktop/ARS/yip_jama_psych/use these/neg_mat_both.csv'

# Path and name of file for results 
output_results_path <- '/Users/bobkohler/Desktop/ARS/ABCD_Projects/alcohol_cpm/'
output_results_filename <- "results_externalmask_Y3NEGEXP_sst.rds"



new_matrices <- readRDS(new_matrices_path)
subject_ids  <- names(new_matrices)
n_subj       <- length(new_matrices)
n_nodes      <- nrow(new_matrices[[1]])


# Makes a 3D array of matrices for each subject. This allows for a more intuitive application of external masks 
x_array <- array(NA_real_, dim = c(n_nodes, n_nodes, n_subj))
for (i in seq_len(n_subj)) {
  mat_i <- new_matrices[[i]]
  stopifnot(
    is.matrix(mat_i),
    nrow(mat_i) == n_nodes,
    ncol(mat_i) == n_nodes
  )
  x_array[ , , i] <- mat_i
}



outcomes_raw <- read.csv(filtered_outcome_path, stringsAsFactors = FALSE)
stopifnot(subject_id_var %in% names(outcomes_raw)) #confirms subject_id variable is found 
stopifnot(outcome_var %in% names(outcomes_raw)) # confirms outcome is found 


# extract only the outcome listed above 
outcomes_df <- outcomes_raw %>%
  mutate(subject_id = as.character(.data[[subject_id_var]]),
         outcome    = .data[[outcome_var]]) %>%
  select(subject_id, outcome)

# separate outcome vector 
outcome_vec <- outcomes_df$outcome

# check to make sure subject_id list and outcome_vec have same IDs
names(outcome_vec) <- outcomes_df$subject_id
stopifnot(all(subject_ids %in% names(outcome_vec)))

y_vec <- outcome_vec[subject_ids]
stopifnot(length(y_vec) == n_subj)


pos_raw <- read.table(pos_raw_path,
                      header = FALSE,
                      stringsAsFactors = FALSE)
neg_raw <- read.table(neg_raw_path,
                      header = FALSE,
                      stringsAsFactors = FALSE)

pos_mat <- as.matrix(pos_raw)
neg_mat <- as.matrix(neg_raw)

# Convert 0/1 notation of masks to TRUE/FALSE masks 
pos_edges <- (pos_mat != 0)
neg_edges <- (neg_mat != 0)

# check dimensions 
stopifnot(
  is.matrix(pos_edges), all(dim(pos_edges) == c(n_nodes, n_nodes)),
  is.matrix(neg_edges), all(dim(neg_edges) == c(n_nodes, n_nodes)))



#~~~~Function for applying external mask and calculating relevant statistics~~~~#
apply_external_mask <- function(x_array, y_vec, pos_edges, neg_edges,
                                outcome_type = c("continuous", "binary")) 
{
  outcome_type <- match.arg(outcome_type)
  
  # Basic dimension checks
  stopifnot(is.array(x_array), length(dim(x_array)) == 3)
  n_nodes <- dim(x_array)[1]
  n_subj  <- dim(x_array)[3]
  stopifnot(length(y_vec) == n_subj)
  stopifnot(is.matrix(pos_edges), all(dim(pos_edges) == c(n_nodes, n_nodes)))
  stopifnot(is.matrix(neg_edges), all(dim(neg_edges) == c(n_nodes, n_nodes)))
  
  # Upper-triangular edges 
  upper_tri <- upper.tri(matrix(0, n_nodes, n_nodes))
  edge_mat  <- apply(x_array, 3, function(m) m[upper_tri])
  mask_p_vec <- pos_edges[upper_tri]
  mask_n_vec <- neg_edges[upper_tri]
  
  # Calculate sum strength
  sum_pos  <- colSums(edge_mat[mask_p_vec, , drop = FALSE])
  sum_neg  <- colSums(edge_mat[mask_n_vec, , drop = FALSE])
  sum_both <- sum_pos - sum_neg
  
  # DF with sum strength and outcome for final analyses 
  df_strength <- tibble(
    subject       = seq_len(n_subj),
    pos_strength  = sum_pos,
    neg_strength  = sum_neg,
    both_strength = sum_both,
    outcome       = y_vec)

  
  # Continuous outcome uses pearson 
  if (outcome_type == "continuous") { 
    cor_pos_val  <- cor.test(df_strength$pos_strength,  df_strength$outcome, method = "spearman")
    cor_neg_val  <- cor.test(df_strength$neg_strength,  df_strength$outcome, method = "spearman")
    cor_both_val <- cor.test(df_strength$both_strength, df_strength$outcome, method = "spearman")
    cor_vals_tbl <- tibble(
      type       = c("pos", "neg", "both"),
      r_value  = c(cor_pos_val$estimate, cor_neg_val$estimate, cor_both_val$estimate),
      p_value = c(cor_pos_val$p.value, cor_neg_val$p.value, cor_both_val$p.value))
    
    glm_pos  <- glm(outcome ~ pos_strength,  data = df_strength, family = Gamma(link = "log"))
    glm_neg  <- glm(outcome ~ neg_strength,  data = df_strength, family = Gamma(link = "log"))
    glm_both <- glm(outcome ~ both_strength, data = df_strength, family = Gamma(link = "log"))
    
    coef_pos  <- summary(glm_pos)$coefficients["pos_strength", ]
    coef_neg  <- summary(glm_neg)$coefficients["neg_strength", ]
    coef_both <- summary(glm_both)$coefficients["both_strength", ]
    
    coefficients_tbl <- tibble(
      term      = c("pos_strength", "neg_strength", "both_strength"),
      estimate  = c(coef_pos["Estimate"],  coef_neg["Estimate"],  coef_both["Estimate"]),
      std_error = c(coef_pos["Std. Error"], coef_neg["Std. Error"], coef_both["Std. Error"]),
      z_value   = c(coef_pos["t value"],    coef_neg["t value"],    coef_both["t value"]),
      p_value   = c(coef_pos["Pr(>|t|)"],   coef_neg["Pr(>|t|)"],   coef_both["Pr(>|t|)"]))
    
    
    return(list(
      strength     = df_strength,
      glm_models   = list(pos = glm_pos, neg = glm_neg, both = glm_both),
      coefficients = coefficients_tbl,
      cor_vals     = cor_vals_tbl
    ))
  }
  
  
  # Binary outcome uses logistic regression & t-test  
  if (outcome_type == "binary") {
    df_strength <- df_strength %>%
      mutate(outcome = factor(as.character(outcome), levels = c("0","1")))
  
  # Logistic 
  glm_pos_binary  <- glm(outcome ~ pos_strength,  data = df_strength, family = "binomial")
  glm_neg_binary  <- glm(outcome ~ neg_strength,  data = df_strength, family = "binomial")
  glm_both_binary <- glm(outcome ~ both_strength, data = df_strength, family = "binomial")
  
  coef_pos_binary  <- summary(glm_pos_binary)$coefficients["pos_strength", ]
  coef_neg_binary <- summary(glm_neg_binary)$coefficients["neg_strength", ]
  coef_both_binary <- summary(glm_both_binary)$coefficients["both_strength", ]
  
  coefficients_glm_binary <- tibble(
    term      = c("pos_strength", "neg_strength", "both_strength"),
    estimate  = c(coef_pos_binary["Estimate"],  coef_neg_binary["Estimate"],  coef_both_binary["Estimate"]),
    std_error = c(coef_pos_binary["Std. Error"], coef_neg_binary["Std. Error"], coef_both_binary["Std. Error"]),
    z_value   = c(coef_pos_binary["z value"],    coef_neg_binary["z value"],    coef_both_binary["z value"]),
    p_value   = c(coef_pos_binary["Pr(>|z|)"],   coef_neg_binary["Pr(>|z|)"],   coef_both_binary["Pr(>|z|)"]))
  
  # AUC 
  roc_pos_binary  <- roc(df_strength$outcome, predict(glm_pos_binary,  type = "response"))
  roc_neg_binary  <- roc(df_strength$outcome, predict(glm_neg_binary,  type = "response"))
  roc_both_binary <- roc(df_strength$outcome, predict(glm_both_binary, type = "response"))
  
  auc_tbl <- tibble(
    type = c("pos", "neg", "both"),
    auc  = c(as.numeric(auc(roc_pos_binary)), as.numeric(auc(roc_neg_binary)), as.numeric(auc(roc_both_binary))))
  
  # t-test
  t_pos  <- t.test(pos_strength  ~ outcome, data = df_strength)
  t_neg  <- t.test(neg_strength  ~ outcome, data = df_strength)
  t_both <- t.test(both_strength ~ outcome, data = df_strength)
  
  coefficients_ttest_binary <- tibble(
    term      = c("pos_strength", "neg_strength", "both_strength"),
    estimate  = c(t_pos["estimate"],  t_neg["estimate"],  t_both["estimate"]),
    std_error = c(t_pos["Std. Error"], t_neg["Std. Error"], t_both["Std. Error"]),
    z_value   = c(t_pos["z value"],    t_neg["z value"],    t_both["z value"]),
    p_value   = c(t_pos["p.value"],   t_neg["p.value"],   t_both["p.value"]))
  
  
  return(list(
    strength   = df_strength,
    glm_models = list(pos = glm_pos_binary, neg = glm_neg_binary, both = glm_both_binary),
    auc        = auc_tbl,
    coefficients_glm = coefficients_glm_binary,
    coefficients_ttest = coefficients_ttest_binary,
    t_tests    = list(pos = t_pos, neg = t_neg, both = t_both)))
  }
}



# Function will run as long as appropriate information is filled in at the top 
results_externalmask <- apply_external_mask(x_array      = x_array, 
                                            y_vec        = y_vec, 
                                            pos_edges    = pos_edges, 
                                            neg_edges    = neg_edges, 
                                            outcome_type = outcome_type)



oldresults <- results_externalmask
# saveRDS(results_externalmask,
#         file.path(output_results_path, output_results_filename)) # This will save the masked matrices thus generating a large file 





