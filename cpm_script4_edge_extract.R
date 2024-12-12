library(R.matlab)
library(tidyverse)

write_path <- '/Users/bobkohler/Desktop/hormone_cpm/hormone_cpm_output/subsample_sex/subsample1_sex_logistic/'
read_path <- '/Users/bobkohler/Desktop/hormone_cpm/hormone_cpm_output/subsample_sex/subsample1_sex_logistic/'

#~~READ IN NODE LIST~~#
lst_rois <- readLines("/Users/bobkohler/Desktop/hormone_cpm/shen268/lst_nodes_orig.txt")

#~~NUMBER OF NODES SHOULD BE 268~~#
num_nodes <- length(lst_rois)
print(num_nodes)

#~~SET THRESHOLD FOR SIGNIFICANT EDGES~~#
thresh <- 1.00

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~POSITIVE EDGE EXTRACTION~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
data_pos <- as.matrix(read.table(file.path(read_path, "network_pos_consensus.txt")))
print(dim(data_pos))

#~~THRESHOLD APPLIED TO POSITIVE NETWORK AND WRITE TABLE~~#
data_thresh_pos <- (data_pos >= thresh) * 1

write.table(
  data_thresh_pos,
  file = file.path(read_path, paste0("net_pos_thresh_", thresh, ".txt")),
  row.names = FALSE,
  col.names = FALSE)


#~~RENAME COLUMNS FOR EASY PLOTTING~~#
pos_threshold_df <- as.data.frame(as.table(data_thresh_pos))
colnames(pos_threshold_df) <- c("Row", "Column", "Value")

#~~GET FREQUENCY MATRIX FOR PLOT~~#
pos_threshold_df <- pos_threshold_df %>%
  mutate(Sum = as.numeric(Row) + as.numeric(Column))

pos_threshold_df$Row <- as.numeric(pos_threshold_df$Row)
pos_threshold_df$Column <- as.numeric(pos_threshold_df$Column)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~POSITIVE EDGE PLOT BY THRESH~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ggplot(pos_threshold_df, 
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
iu_pos <- which(upper.tri(data_thresh_pos), arr.ind = TRUE)

lst_sig_edges_pos <- list()
j_pos <- 0

for (idx in seq_len(nrow(iu_pos))) {
  x_pos <- iu_pos[idx, 1]
  y_pos <- iu_pos[idx, 2]
  j_pos <- j_pos + 1
  if (data_thresh_pos[x_pos, y_pos] == 1) {
    lst_sig_edges_pos <- append(lst_sig_edges_pos, list(c(lst_rois[x_pos], lst_rois[y_pos])))
  }
}

print(j_pos)
print(length(lst_sig_edges_pos))
print(sum(data_thresh_pos) / 2)

#~~SAVE SIGNIFICANT POSITIVE EDGES AFTER APPLYING THRESHOLD~~#
write_conn_pos <- file(file.path(write_path, paste0("sig_edges_pos_thresh", thresh, ".txt")), "w")
for (edge in lst_sig_edges) {
  writeLines(paste(edge, collapse = "\t"), write_conn_pos)
}
close(write_conn_pos)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~NEGATIVE EDGE EXTRACTION~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
data_negative <- as.matrix(read.table(file.path(read_path, "network_neg_consensus.txt")))

#~~THRESHOLD APPLIED TO POSITIVE NETWORK AND WRITE TABLE~~#
data_thresh_neg <- (data_negative >= thresh) * 1

write.table(
  data_thresh_neg,
  file = file.path(read_path, paste0("net_neg_thresh_", thresh, ".txt")),
  row.names = FALSE,
  col.names = FALSE)

#~~RENAME COLUMNS FOR EASY PLOTTING~~#
neg_threshold_df <- as.data.frame(as.table(data_thresh_neg))
colnames(neg_threshold_df) <- c("Row", "Column", "Value")

#~~GET FREQUENCY MATRIX FOR PLOT~~#
neg_threshold_df <- neg_threshold_df %>%
  mutate(Sum = as.numeric(Row) + as.numeric(Column))

neg_threshold_df$Row <- as.numeric(neg_threshold_df$Row)
neg_threshold_df$Column <- as.numeric(neg_threshold_df$Column)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~NEGATIVE EDGE PLOT BY THRESH~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ggplot(neg_threshold_df, 
       aes(x = Column, y = Row, fill = factor(Value))) +
  geom_tile(lwd = 1.5) +
  scale_x_continuous(breaks = seq(0, 268, by = 50)) + # Display all x-axis labels
  scale_y_continuous(breaks = seq(0, 268, by = 50)) + # Display all y-axis labels
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
iu_neg <- which(upper.tri(data_thresh_neg), arr.ind = TRUE)

lst_sig_edges_neg <- list()
j_neg <- 0
for (idx in seq_len(nrow(iu_neg))) {
  x_neg <- iu_neg[idx, 1]
  y_neg <- iu_neg[idx, 2]
  j_neg <- j + 1
  if (data_thresh_neg[x_neg, y_neg] == 1) {
    lst_sig_edges_neg <- append(lst_sig_edges_neg, list(c(lst_rois[x_neg], lst_rois[y_neg])))
  }
}

print(j_neg)
print(length(lst_sig_edges_neg))
print(sum(data_thresh_neg) / 2)

#~~SAVE SIGNIFICANT NEGATIVE EDGES AFTER APPLYING THRESHOLD~~#
write_conn_neg <- file(file.path(write_path, paste0("sig_edges_neg_thresh", thresh, ".txt")), "w")
for (pair in lst_sig_edges_neg) {
  writeLines(paste(pair, collapse = "\t"), write_conn_neg)
}

close(write_conn_neg)


