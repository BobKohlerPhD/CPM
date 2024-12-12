library(jsonlite)
library(dplyr)
library(purrr)
library(parallel)
library(R.matlab)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~UPDATE PATHS TO SCRIPTS WITH REQUIRED CPM FUNCTIONS~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
source("/Users/bobkohler/Desktop/ARS/CPM Scriptscpm_script0_model_functions.R")
source("/Users/bobkohler/Desktop/ARS/CPM Scripts/cpm_script0_run_functions.R")
source("/Users/bobkohler/Desktop/ARS/CPM Scripts/cpm_script0_util_functions.R") 

#~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~CPM PIPELINE FUNCTION~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~#
run_cpm_pipeline <- function(json_file, mat_file_path, mat_file_name, output_dir) {
  
  #~~1. JSON CONFIG LOAD~~#
  config <- fromJSON(json_file)
  
  k <- config$k
  p_thresh <- config$p_thresh
  zscore <- config$zscore
  mode <- config$mode
  num_iterations <- config$num_iter
  
  #~~OUTPUT DIRECTORY~~#
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  #~~2. LOAD .MAT FILE~~#
  data <- read_matlab_mat(mat_file_path, mat_file_name)
  x <- data$x
  y <- data$y
  lst_subj <- data$lst_subjectkey
  
  #~~SAVE SUBJECT ID LIST~~#
  writeLines(lst_subj, file.path(output_dir, "subject_keys.txt"))
  
  #~~3, PARALLELIZE CPM~~#
  cl <- makeCluster(detectCores() - 1)  # Use available cores
  clusterExport(cl, c("run_cpm_thread", "x", "y", "k", "p_thresh", "zscore", "mode", "output_dir"))
  
  #~~RUN CPM IN PARALLEL~~#
  parLapply(cl, seq_len(num_iterations), function(iter) {
    run_cpm_thread(y, iter, x, k, output_dir, p_thresh, zscore, mode)
  })
  stopCluster(cl)
  
  #~~4. COMBINE OUTPUT FILES~~#
  message("Combining results...")
  prediction_files <- list.files(output_dir, pattern = "y_prediction_iter.*\\.csv", full.names = TRUE)
  combined_results <- map_dfr(prediction_files, read_csv)
  
  #~~SAVE OUTPUT~~#
  write_csv(combined_results, file.path(output_dir, "combined_predictions.csv"))
  
  message("CPM pipeline completed. Results saved in ", output_dir)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~SET DIRECTORY PATHS AND RUN CPM~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
json_file <- "path_to_config.json"
mat_file_path <- "path_to_mat_directory"
mat_file_name <- "mat_file_name.mat"
output_dir <- "path_to_output_directory"

#~~~~~~~#
#~~RUN~~#
#~~~~~~~#
run_cpm_pipeline(json_file, mat_file_path, mat_file_name, output_dir)
