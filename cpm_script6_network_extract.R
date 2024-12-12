library(tidyverse)
library(readr)
library(R.matlab)

read_path <- '/Users/bobkohler/Desktop/hormone_cpm/hormone_cpm_output/subsample_sex/subsample1_sex_logistic/'

#~~GET LIST OF MODULE (NETWORK) TEXT FILES IN DIRECTORY~~#
module_files <- list.files(path = file.path(read_path, "canonical_networks"), 
                           pattern = "module.*_comp_realvalued\\.txt", 
                           full.names = TRUE)

module_names <- module_files[!str_detect(module_files, "modules_comp")]
module_names <- module_names[order(as.numeric(str_extract(basename(module_names), "\\d+")))]
module_names


#~~CONFIRM THERE ARE 10~~#
print(length(module_names))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~SET THRESHOLD FOR SIGNIFICANT EDGES~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
thresh <- c(1.00) #USE SAME VALUE AS PREVIOUS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
module_size <- c()

#~~FOR LOOP TO GET MODULE SIZES~~#
for (fn in module_names) {
  lst <- read_lines(fn) 
  module_size <- c(module_size, length(lst)) 
}

#~~PRINT SIZES~~#
print(length(module_size)) #SHOULD BE 10 FOR SHEN
print(sum(module_size)) #SHOULD BE 268 FOR SHEN

#~~READ IN NETWORK TEXT FILES AS LIST TO GET SIG NODES IN EACH NETWORK FROM CPM~~#
lst_rois <- read_lines(file.path(read_path, "canonical_networks/modules_comp_realvalued.txt"))
print(length(lst_rois))

#~~RENAME NETWORKS~~#
lst_nets <- c('medial frontal', 'frontoarietal', 'default', 'motor/sensory', 'visual', 
              'visual b', 'visual association', 'salience', 'subcortical', 'brainstem/cerebellum')
print(length(lst_nets))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~READ IN CANONICAL NETWORK FOR POSITIVE EDGES~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
data_pos <- read.delim(file.path(read_path, "canonical_networks/correlation_canonical_rearr_pos.txt"),
                   header = FALSE, stringsAsFactors = FALSE)

#~~RESTRUCTURE MATRIX AFTER READING IN -- REMOVE ROW ID FROM MATRIX AND MAKE 268*268
rows_pos <- data_pos$V1
processed_pos <- do.call(rbind, lapply(rows_pos, function(rows_pos) {
  values <- strsplit(rows_pos, " ")[[1]]  # split  row by space
  values[-1]  #remove  first value
}))

result_pos <- as.data.frame(processed_pos, stringsAsFactors = FALSE)
data_matrix_pos <- data.frame(lapply(result_pos, as.numeric))


#~~CREATE EMPTY DF WITH NETWORK NAMES~~#
df_pos <- matrix(NA, nrow = 10, ncol = 10, dimnames = list(lst_nets, lst_nets)) %>% as.data.frame()
df_ratio_pos <- matrix(NA, nrow = 10, ncol = 10, dimnames = list(lst_nets, lst_nets)) %>% as.data.frame()


#~~APPLY THRESHOLD TO MATRICES AND IDENTIFY SIGNIFICANT EDGES IN EACH NETWORK~~#
data_thresh_pos <- ifelse(data_matrix_pos >= thresh, 1, 0)

for (i in seq_len(10)) {
  s1 <- module_size[i]
  l1 <- sum(module_size[1:(i-1)])
  h1 <- l1 + s1
  net1 <- lst_nets[i]
  
  for (j in seq_len(10)) {
    s2 <- module_size[j]
    l2 <- sum(module_size[1:(j-1)])
    h2 <- l2 + s2
    net2 <- lst_nets[j]

    num_edges_pos <- sum(data_thresh_pos[(l1+1):h1, (l2+1):h2]) 
    
    if (i == j) {
      num_edges_pos <- num_edges_pos / 2
      network_size <- s1 * (s1 - 1) / 2
    } else {
      network_size <- s1 * s2
    }
    
    df_pos[net1, net2] <- num_edges_pos
    df_ratio_pos[net1, net2] <- num_edges_pos / network_size
  }
}

#~~CONVERT TABLE TO NUMERICS~~#
df_pos <- df_pos %>% mutate(across(everything(), as.numeric))
df_ratio_pos <- df_ratio_pos %>% mutate(across(everything(), as.numeric))

#~~MASK TO CREATE UPPER TRIANGLE FOR EASE OF PLOTTING NETWORKS~~#
mask_pos <- lower.tri(df_pos, diag = FALSE)
df_pos[mask_pos] <- NA

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~READ IN CANONICAL NETWORK FOR NEGATIVE EDGES~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
data_neg <- read.delim(file.path(read_path, "canonical_networks/correlation_canonical_rearr_neg.txt"),
                       header = FALSE, stringsAsFactors = FALSE)

#~~RESTRUCTURE MATRIX AFTER READING IN -- REMOVE ROW ID FROM MATRIX AND MAKE 268*268
rows_neg <- data_neg$V1
processed_neg <- do.call(rbind, lapply(rows_neg, function(rows_neg) {
  values <- strsplit(rows_neg, " ")[[1]] 
}))

result_neg <- as.data.frame(processed_neg, stringsAsFactors = FALSE)
data_matrix_neg <- data.frame(lapply(result_neg, as.numeric))

#~~CREATE EMPTY DF WITH NETWORK NAMES~~#
df_neg <- matrix(NA, nrow = 10, ncol = 10, dimnames = list(lst_nets, lst_nets)) %>% as.data.frame()
df_ratio_neg <- matrix(NA, nrow = 10, ncol = 10, dimnames = list(lst_nets, lst_nets)) %>% as.data.frame()


#~~APPLY THRESHOLD TO MATRICES AND IDENTIFY SIGNIFICANT EDGES IN EACH NETWORK~~#
data_thresh_neg <- ifelse(data_matrix_neg >= thresh, 1, 0)
for (i in seq_len(10)) {
  s1_neg <- module_size[i]
  l1_neg <- sum(module_size[1:(i-1)])
  h1_neg <- l1_neg + s1_neg
  net1_neg <- lst_nets[i]
  
  for (j in seq_len(10)) {
    s2_neg <- module_size[j]
    l2_neg <- sum(module_size[1:(j-1)])
    h2_neg<- l2_neg + s2_neg
    net2_neg <- lst_nets[j]
    
    num_edges_neg <- sum(data_thresh_neg[(l1_neg+1):h1_neg, (l2_neg+1):h2_neg]) 
    
    if (i == j) {
      num_edges_neg <- num_edges_neg / 2
      network_size_neg <- s1_neg * (s1_neg - 1) / 2
    } else {
      network_size_neg <- s1_neg * s2_neg
    }
    
    df_neg[net1_neg, net2_neg] <- num_edges_neg
    df_ratio_neg[net1_neg, net2_neg] <- num_edges_neg / network_size_neg
  }
}

#~~CONVERT TABLE TO NUMERICS~~#
df_neg <- df_neg %>% mutate(across(everything(), as.numeric))
df_ratio_neg <- df_ratio_neg %>% mutate(across(everything(), as.numeric))


#~~MASK TO CREATE UPPER TRIANGLE FOR EASE OF PLOTTING NETWORKS~~#

mask_neg <- lower.tri(df_neg, diag = FALSE)
df_neg[mask_neg] <- NA

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~POSTIVE AND NEGATIVE EDGE PLOT GENERATION~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
x_labels <- c("Medial Frontal", 
              "Frontoparietal", 
              "Default",
              "Motor/Sensory", 
              "Visual",
              "Visual B", 
              "Visual Association", 
              "Salience", 
              "Subcortical", 
              "Brainstem & Cerebellum")

y_labels <- rev(x_labels)

create_heatmap_plot_positive <- function(data, x, y, fill) {
  # Calculate min and max fill values
  min_fill <- min(data[[fill]], na.rm = TRUE)
  max_fill <- max(data[[fill]], na.rm = TRUE)
  
  ggplot(data, aes_string(x = x, y = y, fill = fill)) +
    geom_tile(aes_string(color = sprintf("ifelse(!is.na(%s), 'black', NA)", fill))) + # Dynamically set the border color
    geom_text(aes_string(label = sprintf("ifelse(!is.na(%s), %s, '')", fill, fill),
                         color = sprintf("ifelse(%s > 100, 'white', 'black')", fill)), 
              size = 5.5) +
    scale_fill_gradientn(
      colors = c( "ivory", "#FFE4E1", "lightcoral", "salmon", "indianred2", "brown3", "firebrick4"), # Fine gradient from red to purple
      values = scales::rescale(c(min_fill, 0, max_fill)), # Center gradient around zero
      limits = c(min_fill, max_fill),
      na.value = "white"
    ) + 
    scale_color_identity() +
    scale_x_discrete(labels = x_labels) + 
    scale_y_discrete(labels = y_labels) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(text = element_text(family = "Arial"),
          panel.border = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
          axis.text.y = element_text(angle = 0, vjust = 0.5, face = "bold", size = 12),
          legend.position = "none") +
    coord_fixed()
}
create_heatmap_plot_negative <- function(data, x, y, fill) {
  min_fill <- min(data[[fill]], na.rm = TRUE)
  max_fill <- max(data[[fill]], na.rm = TRUE)
  
  ggplot(data, aes_string(x = x, y = y, fill = fill)) +
    geom_tile(aes_string(color = sprintf("ifelse(!is.na(%s), 'black', NA)", fill))) + # Dynamically set the border color
    geom_text(aes_string(label = sprintf("ifelse(!is.na(%s), %s, '')", fill, fill),
                         color = sprintf("ifelse(%s < 0, 'white', 'black')", fill)),
              size = 5.5) +
    scale_fill_gradientn(
      colors = c("ivory","lightblue", "skyblue2", "dodgerblue", "dodgerblue2", "dodgerblue3", "royalblue3"), # Fine gradient from red to purple
      values = scales::rescale(c(min_fill, 0, max_fill)), # Center gradient around zero
      limits = c(min_fill, max_fill),
      na.value = "white"
    ) +
    scale_color_identity() +
    scale_x_discrete(labels = x_labels) + 
    scale_y_discrete(labels = y_labels) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(text = element_text(family = "Arial"),
          plot.title = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
          axis.text.y = element_text(angle = 0, vjust = 0.5, face = "bold", size = 12),
          legend.position = "none") +
    coord_fixed()
}

#~~GENERATE PLOTS~~#
df_neg_m <- as.matrix(df_neg)
rownames(df_neg) <- colnames(df_neg)
df_neg_melt <- reshape::melt(as.matrix(df_neg), varnames = c("x", "y"), value.name = "value")
df_neg_melt$y <- factor(df_neg_melt$y, levels = rev(colnames(df_neg)))
(negative_edge_plot <- create_heatmap_plot_negative(df_neg_melt, "x", "y", "value"))


df_pos_m <- as.matrix(df_pos)
rownames(df_pos) <- colnames(df_pos)
df_pos_melt <- reshape::melt(as.matrix(df_pos), varnames = c("x", "y"), value.name = "value")
df_pos_melt$y <- factor(df_pos_melt$y, levels = rev(colnames(df_pos)))
(positive_edge_plot <- create_heatmap_plot_positive(df_pos_melt, "x", "y", "value"))
