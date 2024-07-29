# load libraries
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(ggplot2)

# read in data
run11_rawcounts <- read.csv("run11_counts_table.csv", row.names = 1)
run11_rawcounts[] <- round(run11_rawcounts) # DESeq2 requires integer values

run15_rawcounts <- read.csv("run15_counts_table.csv", row.names = 1)
run15_rawcounts[] <- round(run15_rawcounts)

run19_rawcounts <- read.csv("run19_counts_table.csv", row.names = 1)
run19_rawcounts[] <- round(run19_rawcounts)

run23_rawcounts <- read.csv("run23_counts_table.csv", row.names = 1)
run23_rawcounts[] <- round(run23_rawcounts)

run27_rawcounts <- read.csv("run27_counts_table.csv", row.names = 1)
run27_rawcounts[] <- round(run27_rawcounts)

run31_rawcounts <- read.csv("run31_counts_table.csv", row.names = 1)
run31_rawcounts[] <- round(run31_rawcounts)

metadata <- read.delim("DoNoHarm_metadata.tsv")

# reorder rawcounts columns to match metadata
run11_rawcounts <- run11_rawcounts[, match(metadata$ID, colnames(run11_rawcounts))]
run15_rawcounts <- run15_rawcounts[, match(metadata$ID, colnames(run15_rawcounts))]
run19_rawcounts <- run19_rawcounts[, match(metadata$ID, colnames(run19_rawcounts))]
run23_rawcounts <- run23_rawcounts[, match(metadata$ID, colnames(run23_rawcounts))]
run27_rawcounts <- run27_rawcounts[, match(metadata$ID, colnames(run27_rawcounts))]
run31_rawcounts <- run31_rawcounts[, match(metadata$ID, colnames(run31_rawcounts))]

################################################################################
############################### DEFINE FUNCTIONS ###############################
################################################################################

# function to create DESeq2 object and run analysis
run_DESeq2_analysis <- function(counts, meta) {
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = meta,
                                design = ~ timepoint)
  dds <- DESeq(dds)
  results <- results(dds)
  return(as.data.frame(results))
}

# function to generate an MA plot
plot_MA <- function(res, title = "MA Plot") {
  res$color <- with(res, ifelse(padj < 0.05 & log2FoldChange > 1, "green",
                                ifelse(padj < 0.05 & log2FoldChange < -1, 
                                       "red", "black")))
  plot <- ggplot(res, aes(x = log2FoldChange, 
                          y = log10(baseMean), 
                          color = color)) +
    geom_point(alpha = 0.8) +
    scale_color_identity() +
    theme_minimal() +
    labs(x = "Log2 Fold Change", y = "Log10 Mean Counts", title = title) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
    xlim(-10, 10) + ylim(-3, 3) +
    annotate("text", x = 10, y = -3,  
             label = sum(res$color %in% c("green", "red")),  
             hjust = 1, vjust = 1,
             size = 5, color = "black")
  
  # Before using ggsave, make sure 'plot' is a ggplot object and not accidentally a function
  if (!inherits(plot, "ggplot")) stop("plot is not a ggplot object.")
  
  # Save the plot
  plot_filename <- paste0(gsub(" ", "_", title), ".png")
  ggsave(plot_filename, plot, width = 10, height = 8, dpi = 300)
}

process_specific_comparison <- function(comparison_list, counts, metadata) {
  results_list <- list()
  
  for (comparison in comparison_list) {
    # Extract count and metadata for the specific comparison
    conditionA <- paste0("^", comparison[[1]], "_", sep="")
    conditionB <- paste0("^", comparison[[2]], "_", sep="")
    counts_subset <- counts[, grep(paste0(conditionA, "|", conditionB), 
                                   colnames(counts), value = TRUE)]
    meta_subset <- metadata %>% 
      filter(grepl(paste0(conditionA, "|", conditionB), ID))
    
    # Run DESeq2 analysis
    results_df <- run_DESeq2_analysis(counts_subset, meta_subset)
    
    # Plotting
    plot_title <- paste(comparison[[1]], "vs", comparison[[2]])
    plot_MA(results_df, plot_title)
    
    # Append results
    results_list[[plot_title]] <- filter(results_df, padj < 0.05)
  }
  
  return(results_list)
}

# comparisons
comparison_list <- list(
  c("LB_NO_0H", "LB_I_3H"), c("LB_NO_0H", "LB_I_6H"), c("LB_NO_0H", "LB_I_9H"),
  c("LB_NO_0H", "LB_I_24H"), c("LB_NO_0H", "LB_I_48H"), c("LB_I_3H", "LB_I_6H"),
  c("LB_I_3H", "LB_I_9H"), c("LB_I_3H", "LB_I_24H"), c("LB_I_3H", "LB_I_48H"),
  c("LB_I_6H", "LB_I_9H"), c("LB_I_6H", "LB_I_24H"), c("LB_I_6H", "LB_I_48H"),
  c("LB_I_9H", "LB_I_24H"), c("LB_I_9H", "LB_I_48H"), c("LB_I_24H", "LB_I_48H"),
  c("LB_NO_0H", "LB_NI_3H"), c("LB_NO_0H", "LB_NI_6H"), c("LB_NO_0H", "LB_NI_9H"),
  c("LB_NO_0H", "LB_NI_24H"), c("LB_NO_0H", "LB_NI_48H"), c("LB_NI_3H", "LB_NI_6H"),
  c("LB_NI_3H", "LB_NI_9H"), c("LB_NI_3H", "LB_NI_24H"), c("LB_NI_3H", "LB_NI_48H"),
  c("LB_NI_6H", "LB_NI_9H"), c("LB_NI_6H", "LB_NI_24H"), c("LB_NI_6H", "LB_NI_48H"),
  c("LB_NI_9H", "LB_NI_24H"), c("LB_NI_9H", "LB_NI_48H"), c("LB_NI_24H", "LB_NI_48H"),
  c("LB_NO_0H", "DM_I_24H"), c("LB_NO_0H", "DM_I_48H"), c("DM_I_24H", "DM_I_48H"),
  c("LB_NO_0H", "DM_NI_24H"), c("LB_NO_0H", "DM_NI_48H"), c("DM_NI_24H", "DM_NI_48H"))

################################################################################

# time to run
run11_results = process_specific_comparison(comparison_list, run11_rawcounts, metadata)
run15_results = process_specific_comparison(comparison_list, run15_rawcounts, metadata)
run19_results = process_specific_comparison(comparison_list, run19_rawcounts, metadata)
run23_results = process_specific_comparison(comparison_list, run23_rawcounts, metadata)
run27_results = process_specific_comparison(comparison_list, run27_rawcounts, metadata)
run31_results = process_specific_comparison(comparison_list, run31_rawcounts, metadata)
