# Load libraries
library("DESeq2")
library("RColorBrewer")
library("pheatmap")
library("tidyverse")
library("ggplot2")

# count matrix
rawcounts <- read.csv("sequence_counts.csv", row.names = 1)

# read metadata
metadata <- read.delim("DoNoHarm_metadata.tsv")

# check to make sure row names of metadata match column names of rawcounts
all(metadata$ID == colnames(rawcounts))

# create deseq2 object
dds <- DESeqDataSetFromMatrix(countData = rawcounts,
                              colData = metadata,
                              design = ~ timepoint)

# run DESeq2 analysis
dds <- DESeq(dds)
results <- results(dds)
res <- as.data.frame(results)

# plot unfiltered data
res$color <- with(res, ifelse(padj < 0.05 & log2FoldChange > 1, "green",
                              ifelse(padj < 0.05 & log2FoldChange < -1, 
                                     "red", "black")))

volcano_plot <- ggplot(res, aes(x = log2FoldChange, 
                                y = log10(baseMean), 
                                color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "Log10 Mean Counts") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(volcano_plot)

# create data frame of sequences with padj < 0.05 for each comparison
# along with their log2fold change values

################################################################################

# LB, NO, 0H vs LB, I, 3H
LB_NO_0H_vs_LB_I_3H_counts <- rawcounts[, grep("^LB_NO_0H_|^LB_I_3H_", 
                                               colnames(rawcounts), 
                                               value = TRUE)]
LB_NO_0H_vs_LB_I_3H_metadata <- metadata %>% 
  filter(grepl("^LB_NO_0H_|^LB_I_3H_", ID))

LB_NO_0H_vs_LB_I_3H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NO_0H_vs_LB_I_3H_counts,
                         colData = LB_NO_0H_vs_LB_I_3H_metadata,
                         design = ~ timepoint)

LB_NO_0H_vs_LB_I_3H_dds <- DESeq(LB_NO_0H_vs_LB_I_3H_dds)
LB_NO_0H_vs_LB_I_3H_results <- results(LB_NO_0H_vs_LB_I_3H_dds)
LB_NO_0H_vs_LB_I_3H_res <- as.data.frame(LB_NO_0H_vs_LB_I_3H_results)

LB_NO_0H_vs_LB_I_3H_res$color <- 
  with(LB_NO_0H_vs_LB_I_3H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NO_0H_vs_LB_I_3H_volcano_plot <- ggplot(LB_NO_0H_vs_LB_I_3H_res,
                                           aes(x = log2FoldChange,
                                               y = log10(baseMean),
                                               color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NO, 0H vs LB, I, 3H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NO_0H_vs_LB_I_3H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NO_0H_vs_LB_I_3H_volcano_plot)

significant_seqs <- filter(LB_NO_0H_vs_LB_I_3H_res, padj < 0.05)

################################################################################

# LB, NO, 0H vs LB, I, 6H
LB_NO_0H_vs_LB_I_6H_counts <- rawcounts[, grep("^LB_NO_0H_|^LB_I_6H_", 
                                               colnames(rawcounts), 
                                               value = TRUE)]
LB_NO_0H_vs_LB_I_6H_metadata <- metadata %>% 
  filter(grepl("^LB_NO_0H_|^LB_I_6H_", ID))

LB_NO_0H_vs_LB_I_6H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NO_0H_vs_LB_I_6H_counts,
                         colData = LB_NO_0H_vs_LB_I_6H_metadata,
                         design = ~ timepoint)

LB_NO_0H_vs_LB_I_6H_dds <- DESeq(LB_NO_0H_vs_LB_I_6H_dds)
LB_NO_0H_vs_LB_I_6H_results <- results(LB_NO_0H_vs_LB_I_6H_dds)
LB_NO_0H_vs_LB_I_6H_res <- as.data.frame(LB_NO_0H_vs_LB_I_6H_results)

LB_NO_0H_vs_LB_I_6H_res$color <- 
  with(LB_NO_0H_vs_LB_I_6H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NO_0H_vs_LB_I_6H_volcano_plot <- ggplot(LB_NO_0H_vs_LB_I_6H_res,
                                           aes(x = log2FoldChange,
                                               y = log10(baseMean),
                                               color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NO, 0H vs LB, I, 6H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NO_0H_vs_LB_I_6H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NO_0H_vs_LB_I_6H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NO_0H_vs_LB_I_6H_res, padj < 0.05))

################################################################################

# LB, NO, 0H vs LB, I, 9H
LB_NO_0H_vs_LB_I_9H_counts <- rawcounts[, grep("^LB_NO_0H_|^LB_I_9H_", 
                                               colnames(rawcounts), 
                                               value = TRUE)]
LB_NO_0H_vs_LB_I_9H_metadata <- metadata %>% 
  filter(grepl("^LB_NO_0H_|^LB_I_9H_", ID))

LB_NO_0H_vs_LB_I_9H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NO_0H_vs_LB_I_9H_counts,
                         colData = LB_NO_0H_vs_LB_I_9H_metadata,
                         design = ~ timepoint)

LB_NO_0H_vs_LB_I_9H_dds <- DESeq(LB_NO_0H_vs_LB_I_9H_dds)
LB_NO_0H_vs_LB_I_9H_results <- results(LB_NO_0H_vs_LB_I_9H_dds)
LB_NO_0H_vs_LB_I_9H_res <- as.data.frame(LB_NO_0H_vs_LB_I_9H_results)

LB_NO_0H_vs_LB_I_9H_res$color <- 
  with(LB_NO_0H_vs_LB_I_9H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NO_0H_vs_LB_I_9H_volcano_plot <- ggplot(LB_NO_0H_vs_LB_I_9H_res,
                                           aes(x = log2FoldChange,
                                               y = log10(baseMean),
                                               color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NO, 0H vs LB, I, 9H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NO_0H_vs_LB_I_9H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NO_0H_vs_LB_I_9H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NO_0H_vs_LB_I_9H_res, padj < 0.05))

################################################################################

# LB, NO, 0H vs LB, I, 24H
LB_NO_0H_vs_LB_I_24H_counts <- rawcounts[, grep("^LB_NO_0H_|^LB_I_24H_", 
                                               colnames(rawcounts), 
                                               value = TRUE)]
LB_NO_0H_vs_LB_I_24H_metadata <- metadata %>% 
  filter(grepl("^LB_NO_0H_|^LB_I_24H_", ID))

LB_NO_0H_vs_LB_I_24H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NO_0H_vs_LB_I_24H_counts,
                         colData = LB_NO_0H_vs_LB_I_24H_metadata,
                         design = ~ timepoint)

LB_NO_0H_vs_LB_I_24H_dds <- DESeq(LB_NO_0H_vs_LB_I_24H_dds)
LB_NO_0H_vs_LB_I_24H_results <- results(LB_NO_0H_vs_LB_I_24H_dds)
LB_NO_0H_vs_LB_I_24H_res <- as.data.frame(LB_NO_0H_vs_LB_I_24H_results)

LB_NO_0H_vs_LB_I_24H_res$color <- 
  with(LB_NO_0H_vs_LB_I_24H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NO_0H_vs_LB_I_24H_volcano_plot <- ggplot(LB_NO_0H_vs_LB_I_24H_res,
                                           aes(x = log2FoldChange,
                                               y = log10(baseMean),
                                               color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NO, 0H vs LB, I, 24H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NO_0H_vs_LB_I_24H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NO_0H_vs_LB_I_24H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NO_0H_vs_LB_I_24H_res, padj < 0.05))

################################################################################

# LB, NO, 0H vs LB, I, 48H
LB_NO_0H_vs_LB_I_48H_counts <- rawcounts[, grep("^LB_NO_0H_|^LB_I_48H_", 
                                                colnames(rawcounts), 
                                                value = TRUE)]
LB_NO_0H_vs_LB_I_48H_metadata <- metadata %>% 
  filter(grepl("^LB_NO_0H_|^LB_I_48H_", ID))

LB_NO_0H_vs_LB_I_48H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NO_0H_vs_LB_I_48H_counts,
                         colData = LB_NO_0H_vs_LB_I_48H_metadata,
                         design = ~ timepoint)

LB_NO_0H_vs_LB_I_48H_dds <- DESeq(LB_NO_0H_vs_LB_I_48H_dds)
LB_NO_0H_vs_LB_I_48H_results <- results(LB_NO_0H_vs_LB_I_48H_dds)
LB_NO_0H_vs_LB_I_48H_res <- as.data.frame(LB_NO_0H_vs_LB_I_48H_results)

LB_NO_0H_vs_LB_I_48H_res$color <- 
  with(LB_NO_0H_vs_LB_I_48H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NO_0H_vs_LB_I_48H_volcano_plot <- ggplot(LB_NO_0H_vs_LB_I_48H_res,
                                            aes(x = log2FoldChange,
                                                y = log10(baseMean),
                                                color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NO, 0H vs LB, I, 48H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NO_0H_vs_LB_I_48H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NO_0H_vs_LB_I_48H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NO_0H_vs_LB_I_48H_res, padj < 0.05))

################################################################################

# LB, I, 3H vs LB, I, 6H
LB_I_3H_vs_LB_I_6H_counts <- rawcounts[, grep("^LB_I_3H_|^LB_I_6H_", 
                                                colnames(rawcounts), 
                                                value = TRUE)]
LB_I_3H_vs_LB_I_6H_metadata <- metadata %>% 
  filter(grepl("^LB_I_3H_|^LB_I_6H_", ID))

LB_I_3H_vs_LB_I_6H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_I_3H_vs_LB_I_6H_counts,
                         colData = LB_I_3H_vs_LB_I_6H_metadata,
                         design = ~ timepoint)

LB_I_3H_vs_LB_I_6H_dds <- DESeq(LB_I_3H_vs_LB_I_6H_dds)
LB_I_3H_vs_LB_I_6H_results <- results(LB_I_3H_vs_LB_I_6H_dds)
LB_I_3H_vs_LB_I_6H_res <- as.data.frame(LB_I_3H_vs_LB_I_6H_results)

LB_I_3H_vs_LB_I_6H_res$color <- 
  with(LB_I_3H_vs_LB_I_6H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_I_3H_vs_LB_I_6H_volcano_plot <- ggplot(LB_I_3H_vs_LB_I_6H_res,
                                            aes(x = log2FoldChange,
                                                y = log10(baseMean),
                                                color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, I, 3H vs LB, I, 6H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_I_3H_vs_LB_I_6H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_I_3H_vs_LB_I_6H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_I_3H_vs_LB_I_6H_res, padj < 0.05))

################################################################################

# LB, I, 3H vs LB, I, 9H
LB_I_3H_vs_LB_I_9H_counts <- rawcounts[, grep("^LB_I_3H_|^LB_I_9H_", 
                                              colnames(rawcounts), 
                                              value = TRUE)]
LB_I_3H_vs_LB_I_9H_metadata <- metadata %>% 
  filter(grepl("^LB_I_3H_|^LB_I_9H_", ID))

LB_I_3H_vs_LB_I_9H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_I_3H_vs_LB_I_9H_counts,
                         colData = LB_I_3H_vs_LB_I_9H_metadata,
                         design = ~ timepoint)

LB_I_3H_vs_LB_I_9H_dds <- DESeq(LB_I_3H_vs_LB_I_9H_dds)
LB_I_3H_vs_LB_I_9H_results <- results(LB_I_3H_vs_LB_I_9H_dds)
LB_I_3H_vs_LB_I_9H_res <- as.data.frame(LB_I_3H_vs_LB_I_9H_results)

LB_I_3H_vs_LB_I_9H_res$color <- 
  with(LB_I_3H_vs_LB_I_9H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_I_3H_vs_LB_I_9H_volcano_plot <- ggplot(LB_I_3H_vs_LB_I_9H_res,
                                          aes(x = log2FoldChange,
                                              y = log10(baseMean),
                                              color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, I, 3H vs LB, I, 9H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_I_3H_vs_LB_I_9H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_I_3H_vs_LB_I_9H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_I_3H_vs_LB_I_9H_res, padj < 0.05))

################################################################################

# LB, I, 3H vs LB, I, 24H
LB_I_3H_vs_LB_I_24H_counts <- rawcounts[, grep("^LB_I_3H_|^LB_I_24H_", 
                                              colnames(rawcounts), 
                                              value = TRUE)]
LB_I_3H_vs_LB_I_24H_metadata <- metadata %>% 
  filter(grepl("^LB_I_3H_|^LB_I_24H_", ID))

LB_I_3H_vs_LB_I_24H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_I_3H_vs_LB_I_24H_counts,
                         colData = LB_I_3H_vs_LB_I_24H_metadata,
                         design = ~ timepoint)

LB_I_3H_vs_LB_I_24H_dds <- DESeq(LB_I_3H_vs_LB_I_24H_dds)
LB_I_3H_vs_LB_I_24H_results <- results(LB_I_3H_vs_LB_I_24H_dds)
LB_I_3H_vs_LB_I_24H_res <- as.data.frame(LB_I_3H_vs_LB_I_24H_results)

LB_I_3H_vs_LB_I_24H_res$color <- 
  with(LB_I_3H_vs_LB_I_24H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_I_3H_vs_LB_I_24H_volcano_plot <- ggplot(LB_I_3H_vs_LB_I_24H_res,
                                          aes(x = log2FoldChange,
                                              y = log10(baseMean),
                                              color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, I, 3H vs LB, I, 24H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_I_3H_vs_LB_I_24H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_I_3H_vs_LB_I_24H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_I_3H_vs_LB_I_24H_res, padj < 0.05))

################################################################################

# LB, I, 3H vs LB, I, 48H
LB_I_3H_vs_LB_I_48H_counts <- rawcounts[, grep("^LB_I_3H_|^LB_I_48H_", 
                                               colnames(rawcounts), 
                                               value = TRUE)]
LB_I_3H_vs_LB_I_48H_metadata <- metadata %>% 
  filter(grepl("^LB_I_3H_|^LB_I_48H_", ID))

LB_I_3H_vs_LB_I_48H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_I_3H_vs_LB_I_48H_counts,
                         colData = LB_I_3H_vs_LB_I_48H_metadata,
                         design = ~ timepoint)

LB_I_3H_vs_LB_I_48H_dds <- DESeq(LB_I_3H_vs_LB_I_48H_dds)
LB_I_3H_vs_LB_I_48H_results <- results(LB_I_3H_vs_LB_I_48H_dds)
LB_I_3H_vs_LB_I_48H_res <- as.data.frame(LB_I_3H_vs_LB_I_48H_results)

LB_I_3H_vs_LB_I_48H_res$color <- 
  with(LB_I_3H_vs_LB_I_48H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_I_3H_vs_LB_I_48H_volcano_plot <- ggplot(LB_I_3H_vs_LB_I_48H_res,
                                           aes(x = log2FoldChange,
                                               y = log10(baseMean),
                                               color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, I, 3H vs LB, I, 48H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_I_3H_vs_LB_I_48H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_I_3H_vs_LB_I_48H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_I_3H_vs_LB_I_48H_res, padj < 0.05))

################################################################################

# LB, I, 6H vs LB, I, 9H
LB_I_6H_vs_LB_I_9H_counts <- rawcounts[, grep("^LB_I_6H_|^LB_I_9H_", 
                                               colnames(rawcounts), 
                                               value = TRUE)]
LB_I_6H_vs_LB_I_9H_metadata <- metadata %>% 
  filter(grepl("^LB_I_6H_|^LB_I_9H_", ID))

LB_I_6H_vs_LB_I_9H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_I_6H_vs_LB_I_9H_counts,
                         colData = LB_I_6H_vs_LB_I_9H_metadata,
                         design = ~ timepoint)

LB_I_6H_vs_LB_I_9H_dds <- DESeq(LB_I_6H_vs_LB_I_9H_dds)
LB_I_6H_vs_LB_I_9H_results <- results(LB_I_6H_vs_LB_I_9H_dds)
LB_I_6H_vs_LB_I_9H_res <- as.data.frame(LB_I_6H_vs_LB_I_9H_results)

LB_I_6H_vs_LB_I_9H_res$color <- 
  with(LB_I_6H_vs_LB_I_9H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_I_6H_vs_LB_I_9H_volcano_plot <- ggplot(LB_I_6H_vs_LB_I_9H_res,
                                           aes(x = log2FoldChange,
                                               y = log10(baseMean),
                                               color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, I, 6H vs LB, I, 9H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_I_6H_vs_LB_I_9H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_I_6H_vs_LB_I_9H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_I_6H_vs_LB_I_9H_res, padj < 0.05))

################################################################################

# LB, I, 6H vs LB, I, 24H
LB_I_6H_vs_LB_I_24H_counts <- rawcounts[, grep("^LB_I_6H_|^LB_I_24H_", 
                                              colnames(rawcounts), 
                                              value = TRUE)]
LB_I_6H_vs_LB_I_24H_metadata <- metadata %>% 
  filter(grepl("^LB_I_6H_|^LB_I_24H_", ID))

LB_I_6H_vs_LB_I_24H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_I_6H_vs_LB_I_24H_counts,
                         colData = LB_I_6H_vs_LB_I_24H_metadata,
                         design = ~ timepoint)

LB_I_6H_vs_LB_I_24H_dds <- DESeq(LB_I_6H_vs_LB_I_24H_dds)
LB_I_6H_vs_LB_I_24H_results <- results(LB_I_6H_vs_LB_I_24H_dds)
LB_I_6H_vs_LB_I_24H_res <- as.data.frame(LB_I_6H_vs_LB_I_24H_results)

LB_I_6H_vs_LB_I_24H_res$color <- 
  with(LB_I_6H_vs_LB_I_24H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_I_6H_vs_LB_I_24H_volcano_plot <- ggplot(LB_I_6H_vs_LB_I_24H_res,
                                          aes(x = log2FoldChange,
                                              y = log10(baseMean),
                                              color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, I, 6H vs LB, I, 24H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_I_6H_vs_LB_I_24H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_I_6H_vs_LB_I_24H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_I_6H_vs_LB_I_24H_res, padj < 0.05))

################################################################################

# LB, I, 6H vs LB, I, 48H
LB_I_6H_vs_LB_I_48H_counts <- rawcounts[, grep("^LB_I_6H_|^LB_I_48H_", 
                                               colnames(rawcounts), 
                                               value = TRUE)]
LB_I_6H_vs_LB_I_48H_metadata <- metadata %>% 
  filter(grepl("^LB_I_6H_|^LB_I_48H_", ID))

LB_I_6H_vs_LB_I_48H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_I_6H_vs_LB_I_48H_counts,
                         colData = LB_I_6H_vs_LB_I_48H_metadata,
                         design = ~ timepoint)

LB_I_6H_vs_LB_I_48H_dds <- DESeq(LB_I_6H_vs_LB_I_48H_dds)
LB_I_6H_vs_LB_I_48H_results <- results(LB_I_6H_vs_LB_I_48H_dds)
LB_I_6H_vs_LB_I_48H_res <- as.data.frame(LB_I_6H_vs_LB_I_48H_results)

LB_I_6H_vs_LB_I_48H_res$color <- 
  with(LB_I_6H_vs_LB_I_48H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_I_6H_vs_LB_I_48H_volcano_plot <- ggplot(LB_I_6H_vs_LB_I_48H_res,
                                           aes(x = log2FoldChange,
                                               y = log10(baseMean),
                                               color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, I, 6H vs LB, I, 48H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_I_6H_vs_LB_I_48H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_I_6H_vs_LB_I_48H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_I_6H_vs_LB_I_48H_res, padj < 0.05))

################################################################################

# LB, I, 9H vs LB, I, 24H
LB_I_9H_vs_LB_I_24H_counts <- rawcounts[, grep("^LB_I_9H_|^LB_I_24H_", 
                                               colnames(rawcounts), 
                                               value = TRUE)]
LB_I_9H_vs_LB_I_24H_metadata <- metadata %>% 
  filter(grepl("^LB_I_9H_|^LB_I_24H_", ID))

LB_I_9H_vs_LB_I_24H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_I_9H_vs_LB_I_24H_counts,
                         colData = LB_I_9H_vs_LB_I_24H_metadata,
                         design = ~ timepoint)

LB_I_9H_vs_LB_I_24H_dds <- DESeq(LB_I_9H_vs_LB_I_24H_dds)
LB_I_9H_vs_LB_I_24H_results <- results(LB_I_9H_vs_LB_I_24H_dds)
LB_I_9H_vs_LB_I_24H_res <- as.data.frame(LB_I_9H_vs_LB_I_24H_results)

LB_I_9H_vs_LB_I_24H_res$color <- 
  with(LB_I_9H_vs_LB_I_24H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_I_9H_vs_LB_I_24H_volcano_plot <- ggplot(LB_I_9H_vs_LB_I_24H_res,
                                           aes(x = log2FoldChange,
                                               y = log10(baseMean),
                                               color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, I, 9H vs LB, I, 24H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_I_9H_vs_LB_I_24H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_I_9H_vs_LB_I_24H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_I_9H_vs_LB_I_24H_res, padj < 0.05))

################################################################################

# LB, I, 9H vs LB, I, 48H
LB_I_9H_vs_LB_I_48H_counts <- rawcounts[, grep("^LB_I_9H_|^LB_I_48H_", 
                                               colnames(rawcounts), 
                                               value = TRUE)]
LB_I_9H_vs_LB_I_48H_metadata <- metadata %>% 
  filter(grepl("^LB_I_9H_|^LB_I_48H_", ID))

LB_I_9H_vs_LB_I_48H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_I_9H_vs_LB_I_48H_counts,
                         colData = LB_I_9H_vs_LB_I_48H_metadata,
                         design = ~ timepoint)

LB_I_9H_vs_LB_I_48H_dds <- DESeq(LB_I_9H_vs_LB_I_48H_dds)
LB_I_9H_vs_LB_I_48H_results <- results(LB_I_9H_vs_LB_I_48H_dds)
LB_I_9H_vs_LB_I_48H_res <- as.data.frame(LB_I_9H_vs_LB_I_48H_results)

LB_I_9H_vs_LB_I_48H_res$color <- 
  with(LB_I_9H_vs_LB_I_48H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_I_9H_vs_LB_I_48H_volcano_plot <- ggplot(LB_I_9H_vs_LB_I_48H_res,
                                           aes(x = log2FoldChange,
                                               y = log10(baseMean),
                                               color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, I, 9H vs LB, I, 48H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_I_9H_vs_LB_I_48H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_I_9H_vs_LB_I_48H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_I_9H_vs_LB_I_48H_res, padj < 0.05))

################################################################################

# LB, I, 24H vs LB, I, 48H
LB_I_24H_vs_LB_I_48H_counts <- rawcounts[, grep("^LB_I_24H_|^LB_I_48H_", 
                                               colnames(rawcounts), 
                                               value = TRUE)]
LB_I_24H_vs_LB_I_48H_metadata <- metadata %>% 
  filter(grepl("^LB_I_24H_|^LB_I_48H_", ID))

LB_I_24H_vs_LB_I_48H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_I_24H_vs_LB_I_48H_counts,
                         colData = LB_I_24H_vs_LB_I_48H_metadata,
                         design = ~ timepoint)

LB_I_24H_vs_LB_I_48H_dds <- DESeq(LB_I_24H_vs_LB_I_48H_dds)
LB_I_24H_vs_LB_I_48H_results <- results(LB_I_24H_vs_LB_I_48H_dds)
LB_I_24H_vs_LB_I_48H_res <- as.data.frame(LB_I_24H_vs_LB_I_48H_results)

LB_I_24H_vs_LB_I_48H_res$color <- 
  with(LB_I_24H_vs_LB_I_48H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_I_24H_vs_LB_I_48H_volcano_plot <- ggplot(LB_I_24H_vs_LB_I_48H_res,
                                           aes(x = log2FoldChange,
                                               y = log10(baseMean),
                                               color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, I, 24H vs LB, I, 48H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_I_24H_vs_LB_I_48H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_I_24H_vs_LB_I_48H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_I_24H_vs_LB_I_48H_res, padj < 0.05))

################################################################################

# LB, NO, 0H vs LB, NI, 3H
LB_NO_0H_vs_LB_NI_3H_counts <- rawcounts[, grep("^LB_NO_0H_|^LB_NI_3H_", 
                                               colnames(rawcounts), 
                                               value = TRUE)]
LB_NO_0H_vs_LB_NI_3H_metadata <- metadata %>% 
  filter(grepl("^LB_NO_0H_|^LB_NI_3H_", ID))

LB_NO_0H_vs_LB_NI_3H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NO_0H_vs_LB_NI_3H_counts,
                         colData = LB_NO_0H_vs_LB_NI_3H_metadata,
                         design = ~ timepoint)

LB_NO_0H_vs_LB_NI_3H_dds <- DESeq(LB_NO_0H_vs_LB_NI_3H_dds)
LB_NO_0H_vs_LB_NI_3H_results <- results(LB_NO_0H_vs_LB_NI_3H_dds)
LB_NO_0H_vs_LB_NI_3H_res <- as.data.frame(LB_NO_0H_vs_LB_NI_3H_results)

LB_NO_0H_vs_LB_NI_3H_res$color <- 
  with(LB_NO_0H_vs_LB_NI_3H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NO_0H_vs_LB_NI_3H_volcano_plot <- ggplot(LB_NO_0H_vs_LB_NI_3H_res,
                                           aes(x = log2FoldChange,
                                               y = log10(baseMean),
                                               color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NO, 0H vs LB, NI, 3H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NO_0H_vs_LB_NI_3H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NO_0H_vs_LB_NI_3H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NO_0H_vs_LB_NI_3H_res, padj < 0.05))

################################################################################

# LB, NO, 0H vs LB, NI, 6H
LB_NO_0H_vs_LB_NI_6H_counts <- rawcounts[, grep("^LB_NO_0H_|^LB_NI_6H_", 
                                                colnames(rawcounts), 
                                                value = TRUE)]
LB_NO_0H_vs_LB_NI_6H_metadata <- metadata %>% 
  filter(grepl("^LB_NO_0H_|^LB_NI_6H_", ID))

LB_NO_0H_vs_LB_NI_6H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NO_0H_vs_LB_NI_6H_counts,
                         colData = LB_NO_0H_vs_LB_NI_6H_metadata,
                         design = ~ timepoint)

LB_NO_0H_vs_LB_NI_6H_dds <- DESeq(LB_NO_0H_vs_LB_NI_6H_dds)
LB_NO_0H_vs_LB_NI_6H_results <- results(LB_NO_0H_vs_LB_NI_6H_dds)
LB_NO_0H_vs_LB_NI_6H_res <- as.data.frame(LB_NO_0H_vs_LB_NI_6H_results)

LB_NO_0H_vs_LB_NI_6H_res$color <- 
  with(LB_NO_0H_vs_LB_NI_6H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NO_0H_vs_LB_NI_6H_volcano_plot <- ggplot(LB_NO_0H_vs_LB_NI_6H_res,
                                            aes(x = log2FoldChange,
                                                y = log10(baseMean),
                                                color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NO, 0H vs LB, NI, 6H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NO_0H_vs_LB_NI_6H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NO_0H_vs_LB_NI_6H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NO_0H_vs_LB_NI_6H_res, padj < 0.05))

################################################################################

# LB, NO, 0H vs LB, NI, 9H
LB_NO_0H_vs_LB_NI_9H_counts <- rawcounts[, grep("^LB_NO_0H_|^LB_NI_9H_", 
                                                colnames(rawcounts), 
                                                value = TRUE)]
LB_NO_0H_vs_LB_NI_9H_metadata <- metadata %>% 
  filter(grepl("^LB_NO_0H_|^LB_NI_9H_", ID))

LB_NO_0H_vs_LB_NI_9H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NO_0H_vs_LB_NI_9H_counts,
                         colData = LB_NO_0H_vs_LB_NI_9H_metadata,
                         design = ~ timepoint)

LB_NO_0H_vs_LB_NI_9H_dds <- DESeq(LB_NO_0H_vs_LB_NI_9H_dds)
LB_NO_0H_vs_LB_NI_9H_results <- results(LB_NO_0H_vs_LB_NI_9H_dds)
LB_NO_0H_vs_LB_NI_9H_res <- as.data.frame(LB_NO_0H_vs_LB_NI_9H_results)

LB_NO_0H_vs_LB_NI_9H_res$color <- 
  with(LB_NO_0H_vs_LB_NI_9H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NO_0H_vs_LB_NI_9H_volcano_plot <- ggplot(LB_NO_0H_vs_LB_NI_9H_res,
                                            aes(x = log2FoldChange,
                                                y = log10(baseMean),
                                                color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NO, 0H vs LB, NI, 9H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NO_0H_vs_LB_NI_9H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NO_0H_vs_LB_NI_9H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NO_0H_vs_LB_NI_9H_res, padj < 0.05))

################################################################################

# LB, NO, 0H vs LB, NI, 24H
LB_NO_0H_vs_LB_NI_24H_counts <- rawcounts[, grep("^LB_NO_0H_|^LB_NI_24H_", 
                                                colnames(rawcounts), 
                                                value = TRUE)]
LB_NO_0H_vs_LB_NI_24H_metadata <- metadata %>% 
  filter(grepl("^LB_NO_0H_|^LB_NI_24H_", ID))

LB_NO_0H_vs_LB_NI_24H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NO_0H_vs_LB_NI_24H_counts,
                         colData = LB_NO_0H_vs_LB_NI_24H_metadata,
                         design = ~ timepoint)

LB_NO_0H_vs_LB_NI_24H_dds <- DESeq(LB_NO_0H_vs_LB_NI_24H_dds)
LB_NO_0H_vs_LB_NI_24H_results <- results(LB_NO_0H_vs_LB_NI_24H_dds)
LB_NO_0H_vs_LB_NI_24H_res <- as.data.frame(LB_NO_0H_vs_LB_NI_24H_results)

LB_NO_0H_vs_LB_NI_24H_res$color <- 
  with(LB_NO_0H_vs_LB_NI_24H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NO_0H_vs_LB_NI_24H_volcano_plot <- ggplot(LB_NO_0H_vs_LB_NI_24H_res,
                                            aes(x = log2FoldChange,
                                                y = log10(baseMean),
                                                color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NO, 0H vs LB, NI, 24H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NO_0H_vs_LB_NI_24H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NO_0H_vs_LB_NI_24H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NO_0H_vs_LB_NI_24H_res, padj < 0.05))

################################################################################

# LB, NO, 0H vs LB, NI, 48H
LB_NO_0H_vs_LB_NI_48H_counts <- rawcounts[, grep("^LB_NO_0H_|^LB_NI_48H_", 
                                                 colnames(rawcounts), 
                                                 value = TRUE)]
LB_NO_0H_vs_LB_NI_48H_metadata <- metadata %>% 
  filter(grepl("^LB_NO_0H_|^LB_NI_48H_", ID))

LB_NO_0H_vs_LB_NI_48H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NO_0H_vs_LB_NI_48H_counts,
                         colData = LB_NO_0H_vs_LB_NI_48H_metadata,
                         design = ~ timepoint)

LB_NO_0H_vs_LB_NI_48H_dds <- DESeq(LB_NO_0H_vs_LB_NI_48H_dds)
LB_NO_0H_vs_LB_NI_48H_results <- results(LB_NO_0H_vs_LB_NI_48H_dds)
LB_NO_0H_vs_LB_NI_48H_res <- as.data.frame(LB_NO_0H_vs_LB_NI_48H_results)

LB_NO_0H_vs_LB_NI_48H_res$color <- 
  with(LB_NO_0H_vs_LB_NI_48H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NO_0H_vs_LB_NI_48H_volcano_plot <- ggplot(LB_NO_0H_vs_LB_NI_48H_res,
                                             aes(x = log2FoldChange,
                                                 y = log10(baseMean),
                                                 color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NO, 0H vs LB, NI, 48H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NO_0H_vs_LB_NI_48H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NO_0H_vs_LB_NI_48H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NO_0H_vs_LB_NI_48H_res, padj < 0.05))

################################################################################

# LB, NI, 3H vs LB, NI, 6H
LB_NI_3H_vs_LB_NI_6H_counts <- rawcounts[, grep("^LB_NI_3H_|^LB_NI_6H_", 
                                                colnames(rawcounts), 
                                                value = TRUE)]
LB_NI_3H_vs_LB_NI_6H_metadata <- metadata %>% 
  filter(grepl("^LB_NI_3H_|^LB_NI_6H_", ID))

LB_NI_3H_vs_LB_NI_6H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NI_3H_vs_LB_NI_6H_counts,
                         colData = LB_NI_3H_vs_LB_NI_6H_metadata,
                         design = ~ timepoint)

LB_NI_3H_vs_LB_NI_6H_dds <- DESeq(LB_NI_3H_vs_LB_NI_6H_dds)
LB_NI_3H_vs_LB_NI_6H_results <- results(LB_NI_3H_vs_LB_NI_6H_dds)
LB_NI_3H_vs_LB_NI_6H_res <- as.data.frame(LB_NI_3H_vs_LB_NI_6H_results)

LB_NI_3H_vs_LB_NI_6H_res$color <- 
  with(LB_NI_3H_vs_LB_NI_6H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NI_3H_vs_LB_NI_6H_volcano_plot <- ggplot(LB_NI_3H_vs_LB_NI_6H_res,
                                            aes(x = log2FoldChange,
                                                y = log10(baseMean),
                                                color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NI, 3H vs LB, NI, 6H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NI_3H_vs_LB_NI_6H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NI_3H_vs_LB_NI_6H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NI_3H_vs_LB_NI_6H_res, padj < 0.05))

################################################################################

# LB, NI, 3H vs LB, NI, 9H
LB_NI_3H_vs_LB_NI_9H_counts <- rawcounts[, grep("^LB_NI_3H_|^LB_NI_9H_", 
                                                colnames(rawcounts), 
                                                value = TRUE)]
LB_NI_3H_vs_LB_NI_9H_metadata <- metadata %>% 
  filter(grepl("^LB_NI_3H_|^LB_NI_9H_", ID))

LB_NI_3H_vs_LB_NI_9H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NI_3H_vs_LB_NI_9H_counts,
                         colData = LB_NI_3H_vs_LB_NI_9H_metadata,
                         design = ~ timepoint)

LB_NI_3H_vs_LB_NI_9H_dds <- DESeq(LB_NI_3H_vs_LB_NI_9H_dds)
LB_NI_3H_vs_LB_NI_9H_results <- results(LB_NI_3H_vs_LB_NI_9H_dds)
LB_NI_3H_vs_LB_NI_9H_res <- as.data.frame(LB_NI_3H_vs_LB_NI_9H_results)

LB_NI_3H_vs_LB_NI_9H_res$color <- 
  with(LB_NI_3H_vs_LB_NI_9H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NI_3H_vs_LB_NI_9H_volcano_plot <- ggplot(LB_NI_3H_vs_LB_NI_9H_res,
                                            aes(x = log2FoldChange,
                                                y = log10(baseMean),
                                                color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NI, 3H vs LB, NI, 9H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NI_3H_vs_LB_NI_9H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NI_3H_vs_LB_NI_9H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NI_3H_vs_LB_NI_9H_res, padj < 0.05))

################################################################################

# LB, NI, 3H vs LB, NI, 24H
LB_NI_3H_vs_LB_NI_24H_counts <- rawcounts[, grep("^LB_NI_3H_|^LB_NI_24H_", 
                                                colnames(rawcounts), 
                                                value = TRUE)]
LB_NI_3H_vs_LB_NI_24H_metadata <- metadata %>% 
  filter(grepl("^LB_NI_3H_|^LB_NI_24H_", ID))

LB_NI_3H_vs_LB_NI_24H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NI_3H_vs_LB_NI_24H_counts,
                         colData = LB_NI_3H_vs_LB_NI_24H_metadata,
                         design = ~ timepoint)

LB_NI_3H_vs_LB_NI_24H_dds <- DESeq(LB_NI_3H_vs_LB_NI_24H_dds)
LB_NI_3H_vs_LB_NI_24H_results <- results(LB_NI_3H_vs_LB_NI_24H_dds)
LB_NI_3H_vs_LB_NI_24H_res <- as.data.frame(LB_NI_3H_vs_LB_NI_24H_results)

LB_NI_3H_vs_LB_NI_24H_res$color <- 
  with(LB_NI_3H_vs_LB_NI_24H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NI_3H_vs_LB_NI_24H_volcano_plot <- ggplot(LB_NI_3H_vs_LB_NI_24H_res,
                                            aes(x = log2FoldChange,
                                                y = log10(baseMean),
                                                color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NI, 3H vs LB, NI, 24H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NI_3H_vs_LB_NI_24H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NI_3H_vs_LB_NI_24H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NI_3H_vs_LB_NI_24H_res, padj < 0.05))

################################################################################

# LB, NI, 3H vs LB, NI, 48H
LB_NI_3H_vs_LB_NI_48H_counts <- rawcounts[, grep("^LB_NI_3H_|^LB_NI_48H_", 
                                                 colnames(rawcounts), 
                                                 value = TRUE)]
LB_NI_3H_vs_LB_NI_48H_metadata <- metadata %>% 
  filter(grepl("^LB_NI_3H_|^LB_NI_48H_", ID))

LB_NI_3H_vs_LB_NI_48H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NI_3H_vs_LB_NI_48H_counts,
                         colData = LB_NI_3H_vs_LB_NI_48H_metadata,
                         design = ~ timepoint)

LB_NI_3H_vs_LB_NI_48H_dds <- DESeq(LB_NI_3H_vs_LB_NI_48H_dds)
LB_NI_3H_vs_LB_NI_48H_results <- results(LB_NI_3H_vs_LB_NI_48H_dds)
LB_NI_3H_vs_LB_NI_48H_res <- as.data.frame(LB_NI_3H_vs_LB_NI_48H_results)

LB_NI_3H_vs_LB_NI_48H_res$color <- 
  with(LB_NI_3H_vs_LB_NI_48H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NI_3H_vs_LB_NI_48H_volcano_plot <- ggplot(LB_NI_3H_vs_LB_NI_48H_res,
                                             aes(x = log2FoldChange,
                                                 y = log10(baseMean),
                                                 color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NI, 3H vs LB, NI, 48H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NI_3H_vs_LB_NI_48H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NI_3H_vs_LB_NI_48H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NI_3H_vs_LB_NI_48H_res, padj < 0.05))

################################################################################

# LB, NI, 6H vs LB, NI, 9H
LB_NI_6H_vs_LB_NI_9H_counts <- rawcounts[, grep("^LB_NI_6H_|^LB_NI_9H_", 
                                                colnames(rawcounts), 
                                                value = TRUE)]
LB_NI_6H_vs_LB_NI_9H_metadata <- metadata %>% 
  filter(grepl("^LB_NI_6H_|^LB_NI_9H_", ID))

LB_NI_6H_vs_LB_NI_9H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NI_6H_vs_LB_NI_9H_counts,
                         colData = LB_NI_6H_vs_LB_NI_9H_metadata,
                         design = ~ timepoint)

LB_NI_6H_vs_LB_NI_9H_dds <- DESeq(LB_NI_6H_vs_LB_NI_9H_dds)
LB_NI_6H_vs_LB_NI_9H_results <- results(LB_NI_6H_vs_LB_NI_9H_dds)
LB_NI_6H_vs_LB_NI_9H_res <- as.data.frame(LB_NI_6H_vs_LB_NI_9H_results)

LB_NI_6H_vs_LB_NI_9H_res$color <- 
  with(LB_NI_6H_vs_LB_NI_9H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NI_6H_vs_LB_NI_9H_volcano_plot <- ggplot(LB_NI_6H_vs_LB_NI_9H_res,
                                            aes(x = log2FoldChange,
                                                y = log10(baseMean),
                                                color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NI, 6H vs LB, NI, 9H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NI_6H_vs_LB_NI_9H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NI_6H_vs_LB_NI_9H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NI_6H_vs_LB_NI_9H_res, padj < 0.05))

################################################################################

# LB, NI, 6H vs LB, NI, 24H
LB_NI_6H_vs_LB_NI_24H_counts <- rawcounts[, grep("^LB_NI_6H_|^LB_NI_24H_", 
                                                colnames(rawcounts), 
                                                value = TRUE)]
LB_NI_6H_vs_LB_NI_24H_metadata <- metadata %>% 
  filter(grepl("^LB_NI_6H_|^LB_NI_24H_", ID))

LB_NI_6H_vs_LB_NI_24H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NI_6H_vs_LB_NI_24H_counts,
                         colData = LB_NI_6H_vs_LB_NI_24H_metadata,
                         design = ~ timepoint)

LB_NI_6H_vs_LB_NI_24H_dds <- DESeq(LB_NI_6H_vs_LB_NI_24H_dds)
LB_NI_6H_vs_LB_NI_24H_results <- results(LB_NI_6H_vs_LB_NI_24H_dds)
LB_NI_6H_vs_LB_NI_24H_res <- as.data.frame(LB_NI_6H_vs_LB_NI_24H_results)

LB_NI_6H_vs_LB_NI_24H_res$color <- 
  with(LB_NI_6H_vs_LB_NI_24H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NI_6H_vs_LB_NI_24H_volcano_plot <- ggplot(LB_NI_6H_vs_LB_NI_24H_res,
                                            aes(x = log2FoldChange,
                                                y = log10(baseMean),
                                                color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NI, 6H vs LB, NI, 24H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NI_6H_vs_LB_NI_24H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NI_6H_vs_LB_NI_24H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NI_6H_vs_LB_NI_24H_res, padj < 0.05))

################################################################################

# LB, NI, 6H vs LB, NI, 48H
LB_NI_6H_vs_LB_NI_48H_counts <- rawcounts[, grep("^LB_NI_6H_|^LB_NI_48H_", 
                                                 colnames(rawcounts), 
                                                 value = TRUE)]
LB_NI_6H_vs_LB_NI_48H_metadata <- metadata %>% 
  filter(grepl("^LB_NI_6H_|^LB_NI_48H_", ID))

LB_NI_6H_vs_LB_NI_48H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NI_6H_vs_LB_NI_48H_counts,
                         colData = LB_NI_6H_vs_LB_NI_48H_metadata,
                         design = ~ timepoint)

LB_NI_6H_vs_LB_NI_48H_dds <- DESeq(LB_NI_6H_vs_LB_NI_48H_dds)
LB_NI_6H_vs_LB_NI_48H_results <- results(LB_NI_6H_vs_LB_NI_48H_dds)
LB_NI_6H_vs_LB_NI_48H_res <- as.data.frame(LB_NI_6H_vs_LB_NI_48H_results)

LB_NI_6H_vs_LB_NI_48H_res$color <- 
  with(LB_NI_6H_vs_LB_NI_48H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NI_6H_vs_LB_NI_48H_volcano_plot <- ggplot(LB_NI_6H_vs_LB_NI_48H_res,
                                             aes(x = log2FoldChange,
                                                 y = log10(baseMean),
                                                 color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NI, 6H vs LB, NI, 48H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NI_6H_vs_LB_NI_48H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NI_6H_vs_LB_NI_48H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NI_6H_vs_LB_NI_48H_res, padj < 0.05))

################################################################################

# LB, NI, 9H vs LB, NI, 24H
LB_NI_9H_vs_LB_NI_24H_counts <- rawcounts[, grep("^LB_NI_9H_|^LB_NI_24H_", 
                                                 colnames(rawcounts), 
                                                 value = TRUE)]
LB_NI_9H_vs_LB_NI_24H_metadata <- metadata %>% 
  filter(grepl("^LB_NI_9H_|^LB_NI_24H_", ID))

LB_NI_9H_vs_LB_NI_24H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NI_9H_vs_LB_NI_24H_counts,
                         colData = LB_NI_9H_vs_LB_NI_24H_metadata,
                         design = ~ timepoint)

LB_NI_9H_vs_LB_NI_24H_dds <- DESeq(LB_NI_9H_vs_LB_NI_24H_dds)
LB_NI_9H_vs_LB_NI_24H_results <- results(LB_NI_9H_vs_LB_NI_24H_dds)
LB_NI_9H_vs_LB_NI_24H_res <- as.data.frame(LB_NI_9H_vs_LB_NI_24H_results)

LB_NI_9H_vs_LB_NI_24H_res$color <- 
  with(LB_NI_9H_vs_LB_NI_24H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NI_9H_vs_LB_NI_24H_volcano_plot <- ggplot(LB_NI_9H_vs_LB_NI_24H_res,
                                             aes(x = log2FoldChange,
                                                 y = log10(baseMean),
                                                 color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NI, 9H vs LB, NI, 24H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NI_9H_vs_LB_NI_24H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NI_9H_vs_LB_NI_24H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NI_9H_vs_LB_NI_24H_res, padj < 0.05))

################################################################################

# LB, NI, 9H vs LB, NI, 48H
LB_NI_9H_vs_LB_NI_48H_counts <- rawcounts[, grep("^LB_NI_9H_|^LB_NI_48H_", 
                                                 colnames(rawcounts), 
                                                 value = TRUE)]
LB_NI_9H_vs_LB_NI_48H_metadata <- metadata %>% 
  filter(grepl("^LB_NI_9H_|^LB_NI_48H_", ID))

LB_NI_9H_vs_LB_NI_48H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NI_9H_vs_LB_NI_48H_counts,
                         colData = LB_NI_9H_vs_LB_NI_48H_metadata,
                         design = ~ timepoint)

LB_NI_9H_vs_LB_NI_48H_dds <- DESeq(LB_NI_9H_vs_LB_NI_48H_dds)
LB_NI_9H_vs_LB_NI_48H_results <- results(LB_NI_9H_vs_LB_NI_48H_dds)
LB_NI_9H_vs_LB_NI_48H_res <- as.data.frame(LB_NI_9H_vs_LB_NI_48H_results)

LB_NI_9H_vs_LB_NI_48H_res$color <- 
  with(LB_NI_9H_vs_LB_NI_48H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NI_9H_vs_LB_NI_48H_volcano_plot <- ggplot(LB_NI_9H_vs_LB_NI_48H_res,
                                             aes(x = log2FoldChange,
                                                 y = log10(baseMean),
                                                 color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NI, 9H vs LB, NI, 48H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NI_9H_vs_LB_NI_48H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NI_9H_vs_LB_NI_48H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NI_9H_vs_LB_NI_48H_res, padj < 0.05))

################################################################################

# LB, NI, 24H vs LB, NI, 48H
LB_NI_24H_vs_LB_NI_48H_counts <- rawcounts[, grep("^LB_NI_24H_|^LB_NI_48H_", 
                                                 colnames(rawcounts), 
                                                 value = TRUE)]
LB_NI_24H_vs_LB_NI_48H_metadata <- metadata %>% 
  filter(grepl("^LB_NI_24H_|^LB_NI_48H_", ID))

LB_NI_24H_vs_LB_NI_48H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NI_24H_vs_LB_NI_48H_counts,
                         colData = LB_NI_24H_vs_LB_NI_48H_metadata,
                         design = ~ timepoint)

LB_NI_24H_vs_LB_NI_48H_dds <- DESeq(LB_NI_24H_vs_LB_NI_48H_dds)
LB_NI_24H_vs_LB_NI_48H_results <- results(LB_NI_24H_vs_LB_NI_48H_dds)
LB_NI_24H_vs_LB_NI_48H_res <- as.data.frame(LB_NI_24H_vs_LB_NI_48H_results)

LB_NI_24H_vs_LB_NI_48H_res$color <- 
  with(LB_NI_24H_vs_LB_NI_48H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NI_24H_vs_LB_NI_48H_volcano_plot <- ggplot(LB_NI_24H_vs_LB_NI_48H_res,
                                             aes(x = log2FoldChange,
                                                 y = log10(baseMean),
                                                 color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NI, 24H vs LB, NI, 48H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NI_24H_vs_LB_NI_48H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NI_24H_vs_LB_NI_48H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NI_24H_vs_LB_NI_48H_res, padj < 0.05))

################################################################################

# LB, NO, 0H vs DM, I, 24H
LB_NO_0H_vs_DM_I_24H_counts <- rawcounts[, grep("^LB_NO_0H_|^DM_I_24H_", 
                                                 colnames(rawcounts), 
                                                 value = TRUE)]
LB_NO_0H_vs_DM_I_24H_metadata <- metadata %>% 
  filter(grepl("^LB_NO_0H_|^DM_I_24H_", ID))

LB_NO_0H_vs_DM_I_24H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NO_0H_vs_DM_I_24H_counts,
                         colData = LB_NO_0H_vs_DM_I_24H_metadata,
                         design = ~ timepoint)

LB_NO_0H_vs_DM_I_24H_dds <- DESeq(LB_NO_0H_vs_DM_I_24H_dds)
LB_NO_0H_vs_DM_I_24H_results <- results(LB_NO_0H_vs_DM_I_24H_dds)
LB_NO_0H_vs_DM_I_24H_res <- as.data.frame(LB_NO_0H_vs_DM_I_24H_results)

LB_NO_0H_vs_DM_I_24H_res$color <- 
  with(LB_NO_0H_vs_DM_I_24H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NO_0H_vs_DM_I_24H_volcano_plot <- ggplot(LB_NO_0H_vs_DM_I_24H_res,
                                             aes(x = log2FoldChange,
                                                 y = log10(baseMean),
                                                 color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NO, 0H vs DM, I, 24H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NO_0H_vs_DM_I_24H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NO_0H_vs_DM_I_24H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NO_0H_vs_DM_I_24H_res, padj < 0.05))

################################################################################

# LB, NO, 0H vs DM, I, 48H
LB_NO_0H_vs_DM_I_48H_counts <- rawcounts[, grep("^LB_NO_0H_|^DM_I_48H_", 
                                                colnames(rawcounts), 
                                                value = TRUE)]
LB_NO_0H_vs_DM_I_48H_metadata <- metadata %>% 
  filter(grepl("^LB_NO_0H_|^DM_I_48H_", ID))

LB_NO_0H_vs_DM_I_48H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NO_0H_vs_DM_I_48H_counts,
                         colData = LB_NO_0H_vs_DM_I_48H_metadata,
                         design = ~ timepoint)

LB_NO_0H_vs_DM_I_48H_dds <- DESeq(LB_NO_0H_vs_DM_I_48H_dds)
LB_NO_0H_vs_DM_I_48H_results <- results(LB_NO_0H_vs_DM_I_48H_dds)
LB_NO_0H_vs_DM_I_48H_res <- as.data.frame(LB_NO_0H_vs_DM_I_48H_results)

LB_NO_0H_vs_DM_I_48H_res$color <- 
  with(LB_NO_0H_vs_DM_I_48H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NO_0H_vs_DM_I_48H_volcano_plot <- ggplot(LB_NO_0H_vs_DM_I_48H_res,
                                            aes(x = log2FoldChange,
                                                y = log10(baseMean),
                                                color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NO, 0H vs DM, I, 48H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NO_0H_vs_DM_I_48H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NO_0H_vs_DM_I_48H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NO_0H_vs_DM_I_48H_res, padj < 0.05))

################################################################################

# DM, I, 24H vs DM, I, 48H
DM_I_24H_vs_DM_I_48H_counts <- rawcounts[, grep("^DM_I_24H_|^DM_I_48H_", 
                                                colnames(rawcounts), 
                                                value = TRUE)]
DM_I_24H_vs_DM_I_48H_metadata <- metadata %>% 
  filter(grepl("^DM_I_24H_|^DM_I_48H_", ID))

DM_I_24H_vs_DM_I_48H_dds <- 
  DESeqDataSetFromMatrix(countData = DM_I_24H_vs_DM_I_48H_counts,
                         colData = DM_I_24H_vs_DM_I_48H_metadata,
                         design = ~ timepoint)

DM_I_24H_vs_DM_I_48H_dds <- DESeq(DM_I_24H_vs_DM_I_48H_dds)
DM_I_24H_vs_DM_I_48H_results <- results(DM_I_24H_vs_DM_I_48H_dds)
DM_I_24H_vs_DM_I_48H_res <- as.data.frame(DM_I_24H_vs_DM_I_48H_results)

DM_I_24H_vs_DM_I_48H_res$color <- 
  with(DM_I_24H_vs_DM_I_48H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

DM_I_24H_vs_DM_I_48H_volcano_plot <- ggplot(DM_I_24H_vs_DM_I_48H_res,
                                            aes(x = log2FoldChange,
                                                y = log10(baseMean),
                                                color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "DM, I, 24H vs DM, I, 48H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(DM_I_24H_vs_DM_I_48H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(DM_I_24H_vs_DM_I_48H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(DM_I_24H_vs_DM_I_48H_res, padj < 0.05))

################################################################################

# LB, NO, 0H vs DM, NI, 24H
LB_NO_0H_vs_DM_NI_24H_counts <- rawcounts[, grep("^LB_NO_0H_|^DM_NI_24H_", 
                                                 colnames(rawcounts), 
                                                 value = TRUE)]
LB_NO_0H_vs_DM_NI_24H_metadata <- metadata %>% 
  filter(grepl("^LB_NO_0H_|^DM_NI_24H_", ID))

LB_NO_0H_vs_DM_NI_24H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NO_0H_vs_DM_NI_24H_counts,
                         colData = LB_NO_0H_vs_DM_NI_24H_metadata,
                         design = ~ timepoint)

LB_NO_0H_vs_DM_NI_24H_dds <- DESeq(LB_NO_0H_vs_DM_NI_24H_dds)
LB_NO_0H_vs_DM_NI_24H_results <- results(LB_NO_0H_vs_DM_NI_24H_dds)
LB_NO_0H_vs_DM_NI_24H_res <- as.data.frame(LB_NO_0H_vs_DM_NI_24H_results)

LB_NO_0H_vs_DM_NI_24H_res$color <- 
  with(LB_NO_0H_vs_DM_NI_24H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NO_0H_vs_DM_NI_24H_volcano_plot <- ggplot(LB_NO_0H_vs_DM_NI_24H_res,
                                             aes(x = log2FoldChange,
                                                 y = log10(baseMean),
                                                 color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NO, 0H vs DM, NI, 24H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NO_0H_vs_DM_NI_24H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NO_0H_vs_DM_NI_24H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NO_0H_vs_DM_NI_24H_res, padj < 0.05))

################################################################################

# LB, NO, 0H vs DM, NI, 48H
LB_NO_0H_vs_DM_NI_48H_counts <- rawcounts[, grep("^LB_NO_0H_|^DM_NI_48H_", 
                                                 colnames(rawcounts), 
                                                 value = TRUE)]
LB_NO_0H_vs_DM_NI_48H_metadata <- metadata %>% 
  filter(grepl("^LB_NO_0H_|^DM_NI_48H_", ID))

LB_NO_0H_vs_DM_NI_48H_dds <- 
  DESeqDataSetFromMatrix(countData = LB_NO_0H_vs_DM_NI_48H_counts,
                         colData = LB_NO_0H_vs_DM_NI_48H_metadata,
                         design = ~ timepoint)

LB_NO_0H_vs_DM_NI_48H_dds <- DESeq(LB_NO_0H_vs_DM_NI_48H_dds)
LB_NO_0H_vs_DM_NI_48H_results <- results(LB_NO_0H_vs_DM_NI_48H_dds)
LB_NO_0H_vs_DM_NI_48H_res <- as.data.frame(LB_NO_0H_vs_DM_NI_48H_results)

LB_NO_0H_vs_DM_NI_48H_res$color <- 
  with(LB_NO_0H_vs_DM_NI_48H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

LB_NO_0H_vs_DM_NI_48H_volcano_plot <- ggplot(LB_NO_0H_vs_DM_NI_48H_res,
                                             aes(x = log2FoldChange,
                                                 y = log10(baseMean),
                                                 color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "LB, NO, 0H vs DM, NI, 48H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(LB_NO_0H_vs_DM_NI_48H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(LB_NO_0H_vs_DM_NI_48H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(LB_NO_0H_vs_DM_NI_48H_res, padj < 0.05))

################################################################################

# DM, NI, 24H vs DM, NI, 48H
DM_NI_24H_vs_DM_NI_48H_counts <- rawcounts[, grep("^DM_NI_24H_|^DM_NI_48H_", 
                                                 colnames(rawcounts), 
                                                 value = TRUE)]
DM_NI_24H_vs_DM_NI_48H_metadata <- metadata %>% 
  filter(grepl("^DM_NI_24H_|^DM_NI_48H_", ID))

DM_NI_24H_vs_DM_NI_48H_dds <- 
  DESeqDataSetFromMatrix(countData = DM_NI_24H_vs_DM_NI_48H_counts,
                         colData = DM_NI_24H_vs_DM_NI_48H_metadata,
                         design = ~ timepoint)

DM_NI_24H_vs_DM_NI_48H_dds <- DESeq(DM_NI_24H_vs_DM_NI_48H_dds)
DM_NI_24H_vs_DM_NI_48H_results <- results(DM_NI_24H_vs_DM_NI_48H_dds)
DM_NI_24H_vs_DM_NI_48H_res <- as.data.frame(DM_NI_24H_vs_DM_NI_48H_results)

DM_NI_24H_vs_DM_NI_48H_res$color <- 
  with(DM_NI_24H_vs_DM_NI_48H_res, 
       ifelse(padj < 0.05 & log2FoldChange > 1, "green",
              ifelse(padj < 0.05 & log2FoldChange < -1, 
                     "red", "black")))

DM_NI_24H_vs_DM_NI_48H_volcano_plot <- ggplot(DM_NI_24H_vs_DM_NI_48H_res,
                                             aes(x = log2FoldChange,
                                                 y = log10(baseMean),
                                                 color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change",
       y = "Log10 Mean Counts",
       title = "DM, NI, 24H vs DM, NI, 48H") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf,
           label = length(which(DM_NI_24H_vs_DM_NI_48H_res$color == "green")),
           hjust=1, vjust=-1,
           size = 5, color = "green")

print(DM_NI_24H_vs_DM_NI_48H_volcano_plot)

significant_seqs <- rbind(significant_seqs,
                          filter(DM_NI_24H_vs_DM_NI_48H_res, padj < 0.05))
