if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)

setwd("/Users/ioannakoutou/Documents/2nd Semester/Genome Analysis - 1MB462/Project_Files/DEA/")

# Load Count Matrix
countdata <- read.delim("counts.txt", header = TRUE, skip = 1)
rownames(countdata) <- countdata$Geneid
countdata <- countdata[, 7:12]
colnames(countdata) <- c("Control_1", "Control_2", "Control_3",
                         "Heat_treated_42_12h_1", "Heat_treated_42_12h_2", "Heat_treated_42_12h_3")

# Sample MetaData
coldata <- data.frame(
  condition = factor(
    c(rep("control", 3), rep("treated", 3)),
    levels = c("control", "treated")),
  row.names = colnames(countdata)
)

# DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = coldata,
  design = ~ condition
)

# Pre-Filtering 
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)

# Save DESeq2 results
res_df <- as.data.frame(results(dds))
res_df <- res_df[!is.na(res_df$padj), ]
res_df$GeneName <- rownames(res_df)
res_df$Regulation <- ifelse(res_df$padj < 0.05 & res_df$log2FoldChange > 0, "Upregulated",
                            ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < 0, "Downregulated",
                                   "Not Significant"))

write.table(res_df, file = "deseq2_all_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# PCA plot
rld <- rlog(dds, blind = FALSE)
pca_df <- plotPCA(rld, intgroup = "condition", returnData = TRUE) %>%
  rownames_to_column("sample") %>%
  mutate(
    hjust    = if_else(sample == "Control_1" | grepl("Heat_treated", sample), 1, 0),
    nudge_x  = if_else(hjust == 1, -1, 1),
    nudge_y  = if_else(sample == "Heat_treated_42_12h_2", -1, 0)
  )

pca <- ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text(
    aes(label = sample, hjust = hjust),
            nudge_x = pca_df$nudge_x,
            nudge_y = pca_df$nudge_y,
            vjust = 0.5, 
            size = 4) +
  theme_minimal(base_size = 15) +
  labs(
    title    = "PCA of RNA-seq Samples",
    subtitle = "Based on rlog-transformed counts"
  ) +
  theme(
    plot.title    = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
ggsave("PCA_plot.png", plot = pca, width = 8, height = 6, dpi = 100, bg = "white")

# Filtering
res_sig <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 2)
res_sig <- res_sig[order(res_sig$padj), ]
res_sig$Regulation <- ifelse(res_sig$log2FoldChange > 0, 
                             "Upregulated", 
                             "Downregulated")
write.table(res_sig, file = "deseq2_filtered_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Dispersion plot
png("dispersion_plot.png",
    width  = 800,
    height = 600,
    res    = 100,
    bg     = "white")
plotDispEsts(dds)
dev.off()

## Plots
# Bar plot
p_bar <- ggplot(res_sig, aes(x = Regulation, fill = Regulation)) +
  geom_bar(width = 0.5) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Downregulated" = "navy", 
                               "Upregulated" = "firebrick3")) +
  labs(title = "Differentially Expressed Genes",
       x = NULL,
       y = "Number of Genes") +
  ylim(0, 180) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
ggsave("barplot_DEGs.png", p_bar, width = 6, height = 4, dpi = 300, bg = "white")

# MA plot  
p_ma <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = Regulation)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_x_log10() +  
  scale_color_manual(values = c("Upregulated" = "firebrick3", 
                                "Downregulated" = "navy", 
                                "Not Significant" = "lightgrey")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.7) +
  labs(title = "MA Plot of Differential Expression",
       x = "Mean of Normalized Counts (log10)",
       y = "Log2 Fold Change",
       color = "Regulation") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave("MA_plot.png", p_ma, width = 6, height = 5, dpi = 300, bg = "white")

# Volcano plot
res_df$Significance <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 2,
                              ifelse(res_df$log2FoldChange > 0, "Upregulated", "Downregulated"),
                              "Not Significant")

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Upregulated" = "firebrick3", 
                                "Downregulated" = "navy", 
                                "Not Significant" = "lightgrey")) +
  labs(title = "Volcano Plot of Differential Expression",
       subtitle = "Thresholds: |log2FC| > 2, padj < 0.05",
       x = "Log2 Fold Change",
       y = "-Log10(padj)",
       color = "Regulation") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(face = "bold"))
ggsave("volcano_plot.png", p_volcano, width = 6, height = 5, dpi = 300, bg = "white")

# Heatmap
norm_counts <- counts(dds, normalized = TRUE)
sig_genes <- res_sig$GeneName
norm_counts_sig <- norm_counts[sig_genes, ]

colnames(norm_counts_sig) <- c(
  "Control_Rep1", "Control_Rep2", "Control_Rep3",
  "Heat_Treated_Rep1", "Heat_Treated_Rep2", "Heat_Treated_Rep3")

norm_counts_scaled <- t(scale(t(norm_counts_sig)))

annotation_col <- data.frame(
  Condition = factor(
    rep(c("control","treated"), each=3), 
    levels = c("control","treated")))
rownames(annotation_col) <- colnames(norm_counts_scaled)

annotation_colors <- list(
  Condition = c(control = "firebrick3",
                treated = "navy"))
my_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(51)
heatmap_breaks  <- seq(-2,2,length.out=51)

pheatmap(norm_counts_scaled,
         color = my_palette,
         breaks = heatmap_breaks,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         annotation_names_col = FALSE,
         main = "Heatmap of Differentially Expressed Genes",
         fontsize = 10,
         fontsize_row = 8,
         fontsize_col = 8,
         treeheight_row = 30,
         angle_col = 45,
         filename = "heatmap_DEGs.png",
         width = 7,
         height = 8)