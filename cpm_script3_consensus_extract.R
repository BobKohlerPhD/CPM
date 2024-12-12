library(tidyverse)
library(purrr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~EXTRACT EDGE CONCENSUS MAPS: PARAMETER SETTING~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~SET DIRECTORY TO OUTPUT FOLDER~~#
read_path <- '/Users/bobkohler/Desktop/hormone_cpm/hormone_cpm_output/subsample_sex/subsample1_sex_logistic/'

#~~SETUP PARAMETERS~~#
n_repeats=100 #SET BASED ON NUMBER OF CPM REPEATS (USUALLY 100)
k=10 #NUMBER OF FOLDS IN CV 

lst_of_i <- 1:n_repeats
lst_of_k <- 1:k

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~POSITIVE NETWORK PLACEHOLDER~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
network_p <- list()

#~~LOOP OVER ALL ITERATIONS AND FOLDS~~#
for (i in lst_of_i) {
  for (k in lst_of_k) {
    file_path <- sprintf("%s/positive_network_from_training_iter%d_fold_%d.txt", read_path, i, k)
    network_p <- append(network_p, list(as.matrix(read.table(file_path))))
  }
}

#~~LIST TO ARRAY~~#
network_p_array <- array(unlist(network_p), dim = c(nrow(network_p[[1]]), ncol(network_p[[1]]), length(network_p)))

#~~CALCULATE AND SAVE POSITIVE CONSENSUS~~#
network_p_consensus <- apply(network_p_array, c(1, 2), mean)
write.table(network_p_consensus, sprintf("%s/network_pos_consensus.txt", read_path), row.names = FALSE, col.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~NEGATIVE NETWORK PLACEHOLDER~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
network_n <- list()

#~~LOOP OVER ALL ITERATIONS AND FOLDS~~#
for (i in lst_of_i) {
  for (k in lst_of_k) {
    file_path <- sprintf("%s/negative_network_from_training_iter%d_fold_%d.txt", read_path, i, k)
    network_n <- append(network_n, list(as.matrix(read.table(file_path))))
  }
}

#~~LIST TO ARRAY~~#
network_n_array <- array(unlist(network_n), dim = c(nrow(network_n[[1]]), ncol(network_n[[1]]), length(network_n)))

#~~CALCULATE AND SAVE NEGATIVE CONSENSUS~~#
network_n_consensus <- apply(network_n_array, c(1, 2), mean)
write.table(network_n_consensus, sprintf("%s/network_neg_consensus.txt", read_path), row.names = FALSE, col.names = FALSE)
