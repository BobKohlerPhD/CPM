library(tidyverse)

# Paths to subject lists and matrix file 
matrix_sublist_path <- '/Users/bobkohler/Desktop/ARS/ABCD_Projects/abcd_matrices_general/timepoint2/tp2_mid_sublist.csv' #Use if subject list not in matrix file 
outcome_sublist_path    <- '/Users/bobkohler/Desktop/ABCD matrix/tp2_mid_sublist.csv' 
matrices_path           <- '/Users/bobkohler/Desktop/ARS/ABCD_Projects/abcd_matrices_general/timepoint2/tp2_MID_FFD15_matrix.rds'

# Path and RDS file name for filtered matrices 
new_rds_path_name <- '/Users/bobkohler/Desktop/ARS/ABCD_Projects/abcd_matrices_general/timepoint2/tp2_mid_alcohol_matrices.rds'

#Read in matrices and pull IDs 
connectivity_matrices <- read_rds(matrices_path)
goodmatrix_sublist <- data.frame(subject_id = names(connectivity_matrices))

# Read in subject of outcome 
outcome_sublist <- read.csv(outcome_sublist_path)

# Merge connectivity matrices with available outcome data
complete_sublist <- merge(goodmatrix_sublist, outcome_sublist, by = "subject_id")
good_ids <- complete_sublist$subect_id 

# Pull matrices based on 'complete_sublist' 
available_ids <- intersect(good_ids, names(connectivity_matrices)) 
filtered_matrices <- connectivity_matrices[available_ids]

# Check for any IDs in 'complete_sublist' that were not found in the matrix file:
missing_ids <- setdiff(good_ids, available_ids)
if (length(missing_ids) > 0) {
  warning("These IDs were in 'complete_sublist' but not in matrix_file:\n",
          paste(missing_ids, collapse = ", "))
}

# Save filtered matrices 
write_rds(filtered_matrices,new_rds_path_name)
