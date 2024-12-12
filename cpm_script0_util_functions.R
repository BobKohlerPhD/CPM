library(purrr)
library(dplyr)
library(readr)
library(R.matlab)

#~~SUBSAMPLE CPM NO GRID~~#
save_run_outputs_subsample_nogrid <- function(out_path, iter, outputs, y) {
  for (fold in seq_along(outputs$edges_p)) {
    write.table(outputs$edges_p[[fold]], 
                file = file.path(out_path, sprintf("positive_network_from_training_iter%d_fold_%d.txt", iter, fold)), 
                row.names = FALSE, col.names = FALSE)
    write.table(outputs$edges_n[[fold]], 
                file = file.path(out_path, sprintf("negative_network_from_training_iter%d_fold_%d.txt", iter, fold)), 
                row.names = FALSE, col.names = FALSE)
  }
  
  df_y_predict <- tibble(y_pred_both = outputs$y_pred_both, y_actual = y)
  write_csv(df_y_predict, file.path(out_path, sprintf("y_prediction_iter%d.csv", iter)))
  
  df_fit <- tibble(
    both_m = map_dbl(outputs$fit_b, ~ .x$coefficients[2]),
    both_b = map_dbl(outputs$fit_b, ~ .x$coefficients[1])
  )
  write_csv(df_fit, file.path(out_path, sprintf("fit_parameters_iter%d.csv", iter)))
}

#~~SAVE SUBSAMPLE OUTPUT~~#
save_run_outputs_subsample <- function(out_path, iter, outputs, y) {
  save_run_outputs_subsample_nogrid(out_path, iter, outputs, y)  # Save basic outputs
  
  for (fold in seq_along(outputs$best_params)) {
    df_params <- as_tibble(outputs$best_params[[fold]])
    write_csv(df_params, file.path(out_path, sprintf("best_params_iter%d_fold_%d.csv", iter, fold)))
  }
}

#~~IF BEHAVIORAL DATA IS NORMALIZED~~#
y_transform <- function(y, y_norm = "id") {
  if (y_norm == "yj") {
    transformer <- preProcess(y, method = c("YeoJohnson", "center", "scale"))
  } else if (y_norm == "id") {
    transformer <- preProcess(y, method = "none")
  } else if (y_norm == "norm") {
    transformer <- preProcess(y, method = c("center", "scale"))
  } else {
    warning(sprintf("Undefined y_norm %s. Using identity function instead.", y_norm))
    transformer <- preProcess(y, method = "none")
  }
  yn <- predict(transformer, as.data.frame(y))
  list(yn = yn[, 1], transformer = transformer)
}

#~~SAVE AS .MAT~~#
save_matlab_mat <- function(path, matname, x, y, lst_subj) {
  mdict <- list(x = x, y = y, subjectkey = lst_subj)
  writeMat(file.path(path, matname), mdict)
}

#~~READ .MAT~~#
read_matlab_mat <- function(path, matname) {
  mdict <- readMat(file.path(path, matname))
  list(x = mdict$x, y = mdict$y, lst_subjectkey = mdict$subjectkey)
}

#~~SYMMETRY CHECK~~#
check_symmetric <- function(a, rtol = 1e-5, atol = 1e-8) {
  all(abs(a - t(a)) <= (atol + rtol * abs(t(a))), na.rm = TRUE)
}

#~~FILE LIST OF SUBJECT MATRICES~~#
generate_file_list <- function(path, lst_subj, num_roi, num_contrasts, t) {
  sprintf("%s/%s_%dROI_%dcontrasts_corr_matrix_%s.txt", path, lst_subj, num_roi, num_contrasts, t)
}

#~~READ MATRICES AS STACKED ARRAY~~#
read_mats <- function(fn_list) {
  fns <- map(fn_list, ~ read.table(.x, header = FALSE))
  
  if (any(map_lgl(fns, ~ any(is.na(.x))))) {
    stop("ERROR: there are NaNs in the correlation matrices! Please check your data.")
  }
  
  array(unlist(fns), dim = c(nrow(fns[[1]]), ncol(fns[[1]]), length(fns)))
}

#~~ESTIMATOR COEFFICIENT FOR LINEAR/RIDGE/LOGISITIC~~#
return_estimator_coef <- function(est, mode) {
  if (mode %in% c("linear", "logistic")) {
    if (is.null(est)) {
      c(NA, NA)
    } else {
      c(coef(est)[2], coef(est)[1])
    }
  } else if (mode == "ridge") {
    if (is.null(est)) {
      c(NA, NA, NA)
    } else {
      c(coef(est$bestTune)[2], coef(est$bestTune)[1], est$bestTune$alpha)
    }
  } else {
    stop(sprintf("ERROR: mode %s not implemented!", mode))
  }
}

#~~SAVE KFOLD CPM~~#
save_run_outputs <- function(out_path, iter, outputs, y_run, mode = "linear") {
  for (fold in seq_along(outputs$edges_p)) {
    write.table(outputs$edges_p[[fold]], 
                file = file.path(out_path, sprintf("positive_network_from_training_iter%d_fold_%d.txt", iter, fold)), 
                row.names = FALSE, col.names = FALSE)
    write.table(outputs$edges_n[[fold]], 
                file = file.path(out_path, sprintf("negative_network_from_training_iter%d_fold_%d.txt", iter, fold)), 
                row.names = FALSE, col.names = FALSE)
  }
  
  df_y_predict <- tibble(
    y_pred_pos = outputs$y_pred_pos,
    y_pred_neg = outputs$y_pred_neg,
    y_pred_both = outputs$y_pred_both,
    y_actual = y_run
  )
  write_csv(df_y_predict, file.path(out_path, sprintf("y_prediction_iter%d.csv", iter)))
  
  df_fit <- tibble(
    pos_m = map_dbl(outputs$fit_p, ~ return_estimator_coef(.x, mode)[1]),
    pos_b = map_dbl(outputs$fit_p, ~ return_estimator_coef(.x, mode)[2]),
    neg_m = map_dbl(outputs$fit_n, ~ return_estimator_coef(.x, mode)[1]),
    neg_b = map_dbl(outputs$fit_n, ~ return_estimator_coef(.x, mode)[2]),
    both_m = map_dbl(outputs$fit_b, ~ return_estimator_coef(.x, mode)[1]),
    both_b = map_dbl(outputs$fit_b, ~ return_estimator_coef(.x, mode)[2])
  )
  
  if (mode == "ridge") {
    df_fit <- df_fit %>% 
      mutate(
        pos_alpha = map_dbl(outputs$fit_p, ~ return_estimator_coef(.x, mode)[3]),
        neg_alpha = map_dbl(outputs$fit_n, ~ return_estimator_coef(.x, mode)[3]),
        both_alpha = map_dbl(outputs$fit_b, ~ return_estimator_coef(.x, mode)[3])
      )
  }
  
  write_csv(df_fit, file.path(out_path, sprintf("fit_parameters_iter%d.csv", iter)))
}
