library(R.matlab)
library(tidyverse)

read_path1 <- "/Users/bobkohler/Desktop/hormone_cpm/CPM Hormone Output/hormone_cpm_output/subsample_sex/subsample4_sex_logistic/"
read_path2 <- '/Users/bobkohler/Desktop/hormone_cpm/CPM Hormone Output/hormone_cpm_output/subsample_sex/subsample5_sex_logistic/'

write_path1 <- read_path1
write_path2 <- read_path2

#~~NODE LIST~~#
lst_rois <- readLines("/Users/bobkohler/Desktop/hormone_cpm/shen268/lst_nodes_orig.txt")

#~~NUMBER OF NODES SHOULD BE 268~~#
num_nodes <- length(lst_rois)
print(num_nodes)

#~~SET THRESHOLD FOR SIGNIFICANT EDGES DETECTION~~#
thresh <- 1.00

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~EDGE OVERLAP EXTRACT~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~POS AND NEG CONSENSUS MATRIX FOR FIRST CPM~~#
data_1_pos <- as.matrix(read.table(file.path(read_path1, "network_pos_consensus.txt")))
data_1_neg <- as.matrix(read.table(file.path(read_path1, "network_neg_consensus.txt")))


data_1_pos_thresh <- (data_1_pos >= thresh) * 1
data_1_neg_thresh <- (data_1_neg >= thresh) * 1
data_1_both_thresh <- data_1_pos_thresh + data_1_neg_thresh

#~~POS AND NEG CONSENSUS MATRIX FOR SECOND CPM~~#
data_2_pos <- as.matrix(read.table(file.path(read_path2, "network_pos_consensus.txt")))
data_2_neg <- as.matrix(read.table(file.path(read_path2, "network_neg_consensus.txt")))

data_2_pos_thresh <- (data_2_pos >= thresh) * 1
data_2_neg_thresh <- (data_2_neg >= thresh) * 1
data_2_both_thresh <- data_2_pos_thresh + data_2_neg_thresh

#~~OVERLAP DATAFRAMES~~#
data_threshold_overlap_positive <- data_1_pos_thresh * data_2_pos_thresh
data_threshold_overlap_negative <- data_1_neg_thresh * data_2_neg_thresh
data_threshold_overlap_both <- data_1_both_thresh * data_2_both_thresh


#~~RENAME COLUMNS FOR EASY PLOTTING~~#
overlap_threshold_df_pos <- as.data.frame(as.table(data_threshold_overlap_positive))
overlap_threshold_df_neg <- as.data.frame(as.table(data_threshold_overlap_negative))
overlap_threshold_df_both <- as.data.frame(as.table(data_threshold_overlap_both))


colnames(overlap_threshold_df_pos) <- c("Row", "Column", "Value")
colnames(overlap_threshold_df_neg) <- c("Row", "Column", "Value")
colnames(overlap_threshold_df_both) <- c("Row", "Column", "Value")

#~~GET FREQUENCY MATRIX FOR PLOT~~#
overlap_threshold_df_pos <- overlap_threshold_df_pos %>%
  mutate(Sum = as.numeric(Row) + as.numeric(Column))
overlap_threshold_df_pos$Row <- as.numeric(overlap_threshold_df_pos$Row)
overlap_threshold_df_pos$Column <- as.numeric(overlap_threshold_df_pos$Column)


overlap_threshold_df_neg <- overlap_threshold_df_neg %>%
  mutate(Sum = as.numeric(Row) + as.numeric(Column))
overlap_threshold_df_neg$Row <- as.numeric(overlap_threshold_df_neg$Row)
overlap_threshold_df_neg$Column <- as.numeric(overlap_threshold_df_neg$Column)

overlap_threshold_df_both <- overlap_threshold_df_both %>%
  mutate(Sum = as.numeric(Row) + as.numeric(Column))
overlap_threshold_df_both$Row <- as.numeric(overlap_threshold_df_both$Row)
overlap_threshold_df_both$Column <- as.numeric(overlap_threshold_df_both$Column)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~OVERLAPPED EDGES PLOT BY THRESH FOR POSITIVE EDGES~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ggplot(overlap_threshold_df_pos, 
       aes(x = Column, y = Row, fill = factor(Value))) +
  geom_tile(lwd = 1.5) +
  scale_x_continuous(breaks = seq(0, 268, by = 50)) + 
  scale_y_continuous(breaks = seq(0, 268, by = 50)) + 
  coord_cartesian() +
  labs(x = "", y = "") +
  scale_fill_manual(values = c("0" = "white", "1" = "red")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.line = element_line(linewidth = 0),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.75),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.ticks = element_blank())


#~~UPPER TRIANGLE POSITIVE EDGE~~#
iu_overlap_pos <- which(upper.tri(data_threshold_overlap_positive), arr.ind = TRUE)

lst_sig_edges_overlap_pos <- list()
j_overlap_pos <- 0

for (idx in seq_len(nrow(iu_overlap_pos))) {
  x_overlap_pos <- iu_overlap_pos[idx, 1]
  y_overlap_pos <- iu_overlap_pos[idx, 2]
  j_overlap_pos <- j_overlap_pos + 1
  if (data_threshold_overlap_positive[x_overlap_pos, y_overlap_pos] == 1) {
    lst_sig_edges_overlap_pos <- append(lst_sig_edges_overlap_pos, list(c(lst_rois[x_overlap_pos], lst_rois[y_overlap_pos])))
  }
}

print(j_overlap_pos)
print(length(lst_sig_edges_overlap_pos))
print(sum(data_threshold_overlap_positive) / 2)


#~~SAVE SIGNIFICANT POSITIVE EDGES AFTER APPLYING THRESHOLD~~#
# write_conn_overlap <- file(file.path(write_path, paste0("sig_edges_overlap_thresh", thresh, ".txt")), "w")
# for (edge in lst_sig_edges) {
#   writeLines(paste(edge, collapse = "\t"), write_conn_overlap)
# }
# close(write_conn_overlap)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~OVERLAPPED EDGES PLOT BY THRESH FOR NEGATIVE EDGES~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ggplot(overlap_threshold_df_neg, 
       aes(x = Column, y = Row, fill = factor(Value))) +
  geom_tile(lwd = 1.5) +
  scale_x_continuous(breaks = seq(0, 268, by = 50)) + 
  scale_y_continuous(breaks = seq(0, 268, by = 50)) + 
  coord_cartesian() +
  labs(x = "", y = "") +
  scale_fill_manual(values = c("0" = "white", "1" = "blue")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.line = element_line(linewidth = 0),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.75),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.ticks = element_blank())


#~~UPPER TRIANGLE NEGATIVE EDGE~~#
iu_overlap_neg <- which(upper.tri(data_threshold_overlap_negative), arr.ind = TRUE)

lst_sig_edges_overlap_neg <- list()
j_overlap_neg <- 0

for (idx in seq_len(nrow(iu_overlap_neg))) {
  x_overlap_neg <- iu_overlap_neg[idx, 1]
  y_overlap_neg <- iu_overlap_neg[idx, 2]
  j_overlap_neg <- j_overlap_neg + 1
  if (data_threshold_overlap_negative[x_overlap_neg, y_overlap_neg] == 1) {
    lst_sig_edges_overlap_neg <- append(lst_sig_edges_overlap_neg, list(c(lst_rois[x_overlap_neg], lst_rois[y_overlap_neg])))
  }
}

node_df <- do.call(rbind, lapply(lst_sig_edges_overlap_neg, function(x) {
  data.frame(Node1 = x[1], Node2 = x[2], stringsAsFactors = FALSE)
}))

print(j_overlap_neg)
print(length(lst_sig_edges_overlap_neg))
print(sum(data_threshold_overlap_negative) / 2)


#~~SAVE SIGNIFICANT POSITIVE EDGES AFTER APPLYING THRESHOLD~~#
# write_conn_overlap <- file(file.path(write_path, paste0("sig_edges_overlap_thresh", thresh, ".txt")), "w")
# for (edge in lst_sig_edges) {
#   writeLines(paste(edge, collapse = "\t"), write_conn_overlap)
# }
# close(write_conn_overlap)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
x_labels <- c("Medial Frontal", 
              "Frontoparietal", 
              "Default",
              "Motor & Sensory", 
              "Visual",
              "Visual B", 
              "Visual Association", 
              "Salience", 
              "Subcortical", 
              "Brainstem & Cerebellum")

y_labels <- rev(x_labels)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
create_heatmap_plot_positive <- function(data, x, y, fill) {
  min_fill <- min(data[[fill]], na.rm = TRUE)
  max_fill <- max(data[[fill]], na.rm = TRUE)
  ggplot(data, aes_string(x = x, y = y, fill = fill)) +
    geom_tile(aes_string(color = sprintf("ifelse(!is.na(%s), 'black', NA)", fill))) + # Dynamically set the border color
    geom_text(aes_string(label = sprintf("ifelse(!is.na(%s), %s, '')", fill, fill),
                         color = sprintf("ifelse(%s > 100, 'white', 'black')", fill)),  size = 7) +
    scale_fill_gradientn(colors = c( "ivory","orange", "orange2", "darkorange", "chocolate2",  "darkorange3"), # Fine gradient from red to purple
                         values = scales::rescale(c(min_fill, max_fill)), 
                         limits = c(min_fill, max_fill),
                         na.value = "white") + 
    scale_color_identity() +
    scale_x_discrete(labels = x_labels) +
    scale_y_discrete(labels = y_labels) +
    labs(x = "", y = "", title = "Positive Edge Overlap") +
    theme_minimal() +
    theme(text = element_text(family = "Arial"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 16),
          axis.text.y = element_text(angle = 0, vjust = 0.5, face = "bold", size = 16),
          legend.position = "none") +
    coord_fixed()
}

create_heatmap_plot_negative <- function(data, x, y, fill) {
  min_fill <- min(data[[fill]], na.rm = TRUE)
  max_fill <- max(data[[fill]], na.rm = TRUE)
  ggplot(data, aes_string(x = x, y = y, fill = fill)) +
    geom_tile(aes_string(color = sprintf("ifelse(!is.na(%s), 'black', NA)", fill))) + # Dynamically set the border color
    geom_text(aes_string(label = sprintf("ifelse(!is.na(%s), %s, '')", fill, fill),
                         color = sprintf("ifelse(%s < 0, 'white', 'black')", fill)),size =7) +
    scale_fill_gradientn(colors = c("ivory","violet", "darkorchid","purple", "purple1", "purple2"), # Fine gradient from red to purple
                         values = scales::rescale(c(min_fill, max_fill)), # Center gradient around zero
                         limits = c(min_fill, max_fill),
                         na.value = "white" ) +
    scale_color_identity() +
    scale_x_discrete(labels = x_labels) + 
    scale_y_discrete(labels = y_labels) +
    labs(x = "", y = "", title = "Negative Edge Overlap") +
    theme_minimal() +
    theme(text = element_text(family = "Arial"),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 16),
          axis.text.y = element_text(angle = 0, vjust = 0.5, face = "bold", size = 16),
          legend.position = "none") +
    coord_fixed()
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mask_data_1_pos <- read.csv(file.path(read_path1, "mask_pos_thresh_1.csv")) 
mask_data_1_neg <- read.csv(file.path(read_path1, "mask_neg_thresh_1.csv")) 

mask_data_2_pos <- read.csv(file.path(read_path2, "mask_pos_thresh_1.csv")) 
mask_data_2_neg <- read.csv(file.path(read_path2, "mask_neg_thresh_1.csv")) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mask_overlap_pos <- data.frame(x = mask_data_1_pos$x, 
                               y = mask_data_1_pos$y, 
                               value1 = as.numeric(mask_data_1_pos$value),
                               value2 = as.numeric(mask_data_2_pos$value)) %>% 
  mutate(value_abs = abs(value1-value2)) %>%
  mutate(value = value1-value2)

mask_overlap_pos$x <- factor(mask_overlap_pos$x, levels = x_labels)
mask_overlap_pos$y <- factor(mask_overlap_pos$y, levels = rev(x_labels))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
(positive_edge_plot <- create_heatmap_plot_positive(mask_overlap_pos, "x", "y", "value_abs"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mask_overlap_neg <- data.frame(x = mask_data_1_neg$x, 
                               y = mask_data_1_neg$y, 
                               value1 = as.numeric(mask_data_1_neg$value),
                               value2 = as.numeric(mask_data_2_neg$value)) %>% 
  mutate(value_abs = abs(value1-value2)) %>% 
  mutate(value = value1-value2)

mask_overlap_neg$x <- factor(mask_overlap_neg$x, levels = x_labels)
mask_overlap_neg$y <- factor(mask_overlap_neg$y, levels = rev(x_labels))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
(negative_edge_plot <- create_heatmap_plot_negative(mask_overlap_neg, "x", "y", "value_abs"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#















