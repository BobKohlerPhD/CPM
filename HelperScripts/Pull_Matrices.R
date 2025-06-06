library(tidyverse)

# Paths to subject lists and matrix file 
matrix_sublist_path     <- '' #Use if subject list not in matrix file 
outcome_sublist_path    <- '' 
matrices_path           <- ''

# Path and RDS file name for filtered matrices 
filtered_matrix_path    <- ''
filtered_outcome_path   <- ''

filtered_matrix_filename  <- ""
filtered_outcome_filename <- ""



#Read in matrices and pull IDs 
connectivity_matrices <- read_rds(matrices_path)
goodmatrix_sublist <- data.frame(src_subject_id = names(connectivity_matrices)) # make sure subject id variable name matches subject id column in outcome 

# Read in subject of outcome 
outcome_sublist <- read.csv(outcome_sublist_path) %>%  # make sure subject id variable name matches subject id column in matrices 
  select(-X)

# Merge connectivity matrices with available outcome data
complete_sublist <- merge(goodmatrix_sublist, outcome_sublist, by = "src_subject_id")
good_ids <- complete_sublist$src_subject_id 

# Pull matrices based on 'complete_sublist' 
available_ids <- intersect(good_ids, names(connectivity_matrices)) 
filtered_matrices <- connectivity_matrices[available_ids]

# Check for any IDs in 'complete_sublist' that were not found in the matrix file:
missing_ids <- setdiff(good_ids, available_ids)
if (length(missing_ids) > 0) {
  warning("These IDs were in 'complete_sublist' but not in matrix_file:\n",
          paste(missing_ids, collapse = ", "))
}


# Save filtered matrices and outcome file 
write_rds(filtered_matrices,
          path = file.path(filtered_matrix_filename, "tp2_mid_alcohol_matrices.rds"))

write.csv(complete_sublist, 
          path = file.path(filtered_outcome_filename, "tp2_filtered_outcome.csv"))


