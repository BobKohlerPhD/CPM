library(parallel)
library(purrr)
library(jsonlite)
library(stringr)
library(lubridate)
library(fs)

#~~FUNCTION FOR JSON PROCESSING~~#
run_cpm_thread <- function(y_run, i, x, k, out_path, p_thresh, zscore, mode) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~Add your CPM implementation here~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  #~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~Placeholder example:~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~#
  result <- list(y_run = y_run, i = i, summary = summary(x))
  saveRDS(result, file = file.path(out_path, paste0("cpm_result_", i, ".rds")))
}

#~~MAIN CPM FUNCTION~~#
main <- function(json_file, num_jobs) {
  #~~READ JSON~~#
  data <- fromJSON(json_file)
  
  #~~PARAMETER EXTRACT FROM JSON~~#
  today_date <- Sys.Date()
  t <- data$t
  k <- data$k
  p_thresh <- data$p_thresh
  n_repeat <- data$repeat
    num_iter <- data$num_iter
  mat_path <- data$mat_path
  mat_name <- data$mat_name
  zscore <- data$zscore
  mode <- data$mode
  y_norm <- data$y_norm
  base_dir <- data$base_dir
  
  #~~OUTPUT PATH~~#
  out_path <- file.path(base_dir, sprintf("%s_%dfold_p_thresh_%s_repeat%d_iter%d_timepoint_%s_z%d_mode_%s_mat_%s", 
                                          today_date, k, p_thresh, n_repeat, num_iter, t, zscore, mode, 
                                          str_remove(mat_name, "\\.mat")))
  dir_create(out_path)
  
  #~~SAVE RUN SETTINGS FILE~~#
  run_settings <- sprintf("Run Date: %s\nJson file: %s\nTime Point: %s\nNumber of folds k: %d\nP threshold: %s\nNumber of repeats: %d\nNumber of iterations: %d\nPath to mat: %s\nmat name: %s\nz-score training edges: %d\nmode: %s\ny norm method: %s\nOutput path: %s", 
                          today_date, json_file, t, k, p_thresh, n_repeat, num_iter, mat_path, mat_name, zscore, mode, y_norm, out_path)
  writeLines(run_settings, file.path(out_path, "run_settings.txt"))
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~READ IN DATA--UPDATE WITH CORRECT MAT FILE NAME)~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  mat <- readRDS(file.path(mat_path, mat_name)) #~~USE CORRECT READER FOR MAT FILES~~#
  x <- mat$x
  y <- mat$y
  lst_subjectkey <- mat$subjectkey
  
  #~~SAVE SUBJECT ID LIST~~#
  writeLines(lst_subjectkey, file.path(out_path, "lst_subjkey_analyzed.txt"))
  
  #~~ITERATION LIST GENERATION~~#
  lst_of_i <- seq_len(n_repeat + num_iter)
  lst_of_yrun <- map(lst_of_i, function(i) {
    if (i <= n_repeat) {
      y  #~~TRUE VALUES FOR BEHAVIOR~~#
    } else {
      sample(y)  #~~PERMUTED Y VALUES~~#
    }
  })
  
  #~~CPU PROCESSESS~~#
  num_cpu <- detectCores()
  use_cpu <- floor(num_cpu * 0.5)
  if (num_jobs >= use_cpu) {
    warning(sprintf("Using more than half of total CPU (%d). Changing number of jobs to half of CPU (%d).", num_jobs, use_cpu))
    num_proc <- use_cpu
  } else {
    num_proc <- num_jobs
  }
  message(sprintf("Using %d jobs", num_proc))
  
  #~~PARALLELIZE~~#
  cl <- makeCluster(num_proc)
  clusterExport(cl, c("run_cpm_thread", "x", "k", "out_path", "p_thresh", "zscore", "mode"))
  parLapply(cl, seq_along(lst_of_yrun), function(i) {
    run_cpm_thread(lst_of_yrun[[i]], i, x, k, out_path, p_thresh, zscore, mode)
  })
  stopCluster(cl)
}

#~~MAIN CPM FUNCTION - RUN AFTER FUNCTIONS HAVE BEEN UPDATED AS NECESSARY AND CORRECT JSON INFO IS SUPPLIE~~#
args <- commandArgs(trailingOnly = TRUE)
json_file <- args[1]
num_jobs <- as.numeric(args[2])
main(json_file, num_jobs)
