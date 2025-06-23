library(tidyverse)
rds_file <- ""

connectivity_mat <- readRDS(rds_file)
ids       <- names(connectivity_mat)
n_subj    <- length(connectivity_mat)

num_missing_edges <- sapply(connectivity_mat, function(mat) {
  sum(mat == 0, na.rm = TRUE)
})

df_missing <- tibble::tibble(
  subj_id       = ids,
  missing_edges = num_missing_edges)

cat(sum(num_missing_edges > 0),
    "subjects have at least one missing edge.\n")


x_array <- simplify2array(connectivity_mat)
numSubs_missingEdges <- apply(x_array == 0, c(1,2), sum)
edge_df <- as.data.frame(numSubs_missingEdges)

