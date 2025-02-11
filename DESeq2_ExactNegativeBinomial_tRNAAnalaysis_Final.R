# Author: Jerrin Thomas George
# Affilliated Lab: Wiedenheft Lab

### Laod packages

# Following packages will be required through out this Rscript. Install necessary packages

# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages("readr")
# install.packages("ggplot2")
# install.packages("reshape2")
# install.packages("writexl")

library(DESeq2)
library(readr)
library(ggplot2)
library(reshape2)
library(writexl)

# Function to read data
Get_data <- function(path_name, filename) {
  data <- read_tsv(path_name)
  row_num <- nrow(data) - 1
  filename <- data[2:row_num, ] # Extract all data from second row to row_num
  return(filename)
}

# Datasets provided within the corresponding GitHub repository

# Load datasets
Cont1 <- Get_data("Set2-Con1-untreated.abundances.tsv", Cont1)
Cont2 <- Get_data("Set2-Con2-untreated.abundances.tsv", Cont2)
Cont3 <- Get_data("Set2-Con3-untreated.abundances.tsv", Cont3)
Exp1 <- Get_data("Set2-Exp1-untreated.abundances.tsv", Exp1)
Exp2 <- Get_data("Set2-Exp2-untreated.abundances.tsv", Exp2)
Exp3 <- Get_data("Set2-Exp3-untreated.abundances.tsv", Exp3)

# Combine Control and Experiment data
CombinedCont <- merge(Cont1[, c(1, 3)], Cont2[, c(1, 3)], by = "Reference sequence", all = FALSE)
CombinedCont <- merge(CombinedCont, Cont3[, c(1, 3)], by = "Reference sequence", all = FALSE)
colnames(CombinedCont) <- c("Reference sequence", "Control1", "Control2", "Control3")

CombinedExp <- merge(Exp1[, c(1, 3)], Exp2[, c(1, 3)], by = "Reference sequence", all = FALSE)
CombinedExp <- merge(CombinedExp, Exp3[, c(1, 3)], by = "Reference sequence", all = FALSE)
colnames(CombinedExp) <- c("Reference sequence", "Experiment1", "Experiment2", "Experiment3")

Combined_Exp_Cont_All <- merge(CombinedCont[, 1:4], CombinedExp[, 1:4], by = "Reference sequence", all = FALSE)
counts <- Combined_Exp_Cont_All[, 2:7]
rownames(counts) <- Combined_Exp_Cont_All[, 1]

# Create Sample Information
sample_info <- data.frame(
  row.names = colnames(counts),
  condition = factor(c("Control", "Control", "Control", "Experiment", "Experiment", "Experiment"))
)

# Create DESeqDataSet object and run analysis
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_info, design = ~ condition)
dds <- DESeq(dds)

# Extract results
res <- results(dds, alpha = 0.05)
resOrdered <- res[order(res$pvalue), ]
res_df <- as.data.frame(resOrdered)
res_df$padj[is.na(res_df$padj)] <- 1
res_df$Significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.585, "Yes", "No")
write.csv(res_df, file = "DESeq2_tRNA_DE_Results.csv", row.names = TRUE)

# Normalize counts
normalized_counts <- counts(dds, normalized = TRUE)
normalized_Control <- normalized_counts[, sample_info$condition == "Control"]
normalized_Experiment <- normalized_counts[, sample_info$condition == "Experiment"]

# Calculate means and standard deviations
mean_Ctrl <- rowMeans(normalized_Control)
mean_Exp <- rowMeans(normalized_Experiment)
sd_Control <- apply(normalized_Control, 1, sd)
sd_Experiment <- apply(normalized_Experiment, 1, sd)

# Combine into a statistics dataframe
Stats_tRNA <- data.frame(
  tRNA = rownames(normalized_counts),
  Control_Mean = mean_Ctrl,
  Experiment_Mean = mean_Exp,
  Control_SD = sd_Control,
  Experiment_SD = sd_Experiment
)

# Save statistics to Excel
write_xlsx(Stats_tRNA, path = "DESeq2_tRNA_Stats.xlsx")

# Volcano Plot
FC_threshold <- log2(1.5) # Fold-change threshold
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = Significant), size = 3, alpha = 1) +
  scale_color_manual(values = c("black", "red"), labels = c("Not Significant", "Significant")) +
  theme_minimal() +
  ylim(-1, 33) +
  xlim(-1.4, 1) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  ggtitle("Volcano Plot of tRNA Differential Expression using DESeq2") +
  geom_vline(xintercept = c(-FC_threshold, FC_threshold), linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", alpha = 0.5) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.25, "cm"),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(size = 20),
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5)
  )

# Prepare data for bar graph
mean_long <- melt(Stats_tRNA, id.vars = "tRNA", 
                  measure.vars = c("Control_Mean", "Experiment_Mean"),
                  variable.name = "Group",
                  value.name = "Mean_CPM")

sd_long <- melt(Stats_tRNA, id.vars = "tRNA",
                measure.vars = c("Control_SD", "Experiment_SD"),
                variable.name = "Group_SD",
                value.name = "SD_CPM")

mean_long$Group <- ifelse(grepl("Control", mean_long$Group), "Control", "Experiment")
sd_long$Group <- ifelse(grepl("Control", sd_long$Group_SD), "Control", "Experiment")
plot_data <- merge(mean_long, sd_long[, c("tRNA", "Group", "SD_CPM")], by = c("tRNA", "Group"))

# Plot bar graph with error bars
ggplot(plot_data, aes(x = tRNA, y = Mean_CPM, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = Mean_CPM - SD_CPM, ymax = Mean_CPM + SD_CPM),
                width = 0.2, 
                position = position_dodge(width = 0.9)) +
  theme_minimal() +
  xlab("tRNA") +
  ylab("Mean Normalized Counts") +
  ggtitle("Mean Normalized Counts of tRNAs in Control and Experiment Groups (DESeq2)") +
  scale_fill_manual(values = c("Control" = "lightblue", "Experiment" = "salmon")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
