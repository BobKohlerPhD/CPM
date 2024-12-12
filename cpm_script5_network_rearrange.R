library(tidyverse)
library(readr)
library(R.matlab)

base_dir <- '/Users/bobkohler/Desktop/hormone_cpm/hormone_cpm_output/subsample_sex/subsample1_sex_logistic/'
results_dir <- file.path(base_dir, "canonical_networks")

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~READ IN SHEN ATLAS MAPPINGS AND UNLIST .MAT~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
lst_rois <- readMat("/Users/bobkohler/Desktop/hormone_cpm/shen268/lst_rois.mat")

nets <- readMat("/Users/bobkohler/Desktop/hormone_cpm/shen268/nets.mat")
nets <- unlist(nets)

lst_nodes_txt <- readLines("/Users/bobkohler/Desktop/hormone_cpm/shen268/lst_nodes_orig.txt")

nodes<- readMat("/Users/bobkohler/Desktop/hormone_cpm/shen268/nodes.mat")
nodes <- unlist(nodes)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~SET MODEL AND ATLAS PARAMETERS~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
n <- 268 #NODE NUMBER
rep <- 100 #NUMBER OF REPTITIONS FOR CONSENSUS CLUSTERING

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~READ IN POSTIVE AND NEGATIVE NETWORK CONSENSUS~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
network_pos_consensus <- as.matrix(read.delim(file.path(base_dir, "network_pos_consensus.txt"), header = FALSE))
network_neg_consensus <- as.matrix(read.delim(file.path(base_dir, "network_neg_consensus.txt"), header = FALSE))

#~~REARRANGE MAT (POSITIVE)~~#
fn_corr_rearr_mat_pos <- file.path(results_dir, "correlation_canonical_rearr_pos.txt")
correlation_rearr_pos <- network_pos_consensus[, 1]
write.table(correlation_rearr_pos, fn_corr_rearr_mat_pos, row.names = T, col.names = F)

#~~REARRANGE MAT (NEGATIVE)~~#
fn_corr_rearr_mat_neg <- file.path(results_dir, "correlation_canonical_rearr_neg.txt")
correlation_rearr_neg <- network_neg_consensus[,1]
write.table(correlation_rearr_neg, fn_corr_rearr_mat_neg, row.names = FALSE, col.names = FALSE)

#~~REARRANGE ROIS BY NETWORK~~#
lstrois_rearr <- lst_nodes_txt[nodes]

#~~GET NETWORK SIZES~~#
num_module <- max(nets)
module_sizes <- sapply(1:num_module, function(i) sum(nets == i))
module_sizes <- as.character(module_sizes)
modules <- list()

for (i in 1:num_module) {
  tmp_module <- lstrois_rearr[which(nets == i)]
  modules <- c(modules, list(tmp_module))
  out_file <- file.path(results_dir, paste0("module", i, "_comp_realvalued.txt"))
  
 write.table(as.character(tmp_module), file = out_file, row.names = FALSE, col.names = FALSE, quote = FALSE) #WRITES TABLE TO OUTPUT FOLDER
}

