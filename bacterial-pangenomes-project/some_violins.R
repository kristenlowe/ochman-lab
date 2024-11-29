library(tidyverse)
library(gridExtra)
library(ggsignif)

#### M TUBERCULOSIS ####

df <- read_csv("/stor/scratch/Ochman/kristen/pangenome/all_bacterial_species/Mycobacterium_tuberculosis/rep_seq_properties/Mycobacterium_tuberculosis_seq_properties.csv", 
               col_names = FALSE)
colnames(df) <- c("gene_id", "feature", "value")

features <- c("length", "GC", "GC_3rd", "polarAA", "hydrophobicAA", "metabol", 
              "cai", "instability", "disorder", "aa_comp_bias", "transmembrane", "tm_coverage")

wide_df <- df %>% 
  pivot_wider(names_from = feature, values_from = value) %>% # 16756 rows
  dplyr::filter((conservation_percentage >= 0.95 & conservation_percentage <= 1) | conservation_percentage <= 0.05) %>% # 14043 rows
  mutate(conservation_classification = ifelse(conservation_percentage >= 0.95, "core", "cloud")) %>% 
  mutate(transmembrane = ifelse(transmembrane == 0, NA, transmembrane)) %>% 
  mutate(tm_coverage = ifelse(tm_coverage == 0, NA, tm_coverage)) %>% 
  mutate(across(all_of(features), 
                ~ ifelse(. >= quantile(., 0.05, na.rm = TRUE) & . <= quantile(., 0.95, na.rm = TRUE), ., NA)))



# Map features to their descriptive titles
feature_titles <- c(
  "length" = "Length (bp)",
  "GC" = "GC Content (%)",
  "GC_3rd" = "GC 3rd Content (%)",
  "polarAA" = "Polar AA Content (%)",
  "hydrophobicAA" = "Hydrophobic AA Content (%)",
  "metabol" = "Metabolic Cost (ATP/mol AA)",
  "cai" = "Codon Adaptation Index",
  "instability" = "Instability Index",
  "disorder" = "Intrinsic Structural Disorder",
  "aa_comp_bias" = "Amino Acid Composition Bias",
  "transmembrane" = "Transmembrane Domain Count",
  "tm_coverage" = "Transmembrane Coverage"
)

# Updated plot function
plots <- lapply(features, function(feature) {
  ggplot(wide_df, aes(x = conservation_classification, y = .data[[feature]])) +
    geom_violin(aes(fill = conservation_classification), alpha = 1, trim = TRUE) +
    scale_color_manual(values = c("tan1", "skyblue")) +
    scale_fill_manual(values = c("tan1", "skyblue")) +
    geom_signif(
      comparisons = list(c("cloud", "core")),
      map_signif_level = TRUE
    ) +
    theme_minimal() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(x = "", y = feature_titles[feature]) +
    ggtitle(feature_titles[feature])
})

# Save the arranged grid with updated titles
ggsave("Mycobacterium_tuberculosis_violin.png", 
       plot = do.call(grid.arrange, 
                      c(plots, 
                        ncol = 6, 
                        top = "Mycobacterium tuberculosis")), 
       width = 20, height = 14)


#### E COLI ####

df <- read_csv("/stor/scratch/Ochman/kristen/pangenome/all_bacterial_species/Escherichia_coli/rep_seq_properties/Escherichia_coli_seq_properties.csv", 
               col_names = FALSE)
colnames(df) <- c("gene_id", "feature", "value")

wide_df <- df %>% 
  pivot_wider(names_from = feature, values_from = value) %>% # 16756 rows
  dplyr::filter((conservation_percentage >= 0.95 & conservation_percentage <= 1) | conservation_percentage <= 0.05) %>% # 14043 rows
  mutate(conservation_classification = ifelse(conservation_percentage >= 0.95, "core", "cloud")) %>% 
  mutate(transmembrane = ifelse(transmembrane == 0, NA, transmembrane)) %>% 
  mutate(tm_coverage = ifelse(tm_coverage == 0, NA, tm_coverage)) %>% 
  mutate(across(all_of(features), 
                ~ ifelse(. >= quantile(., 0.05, na.rm = TRUE) & . <= quantile(., 0.95, na.rm = TRUE), ., NA)))

plots <- lapply(features, function(feature) {
  ggplot(wide_df, aes(x = conservation_classification, y = .data[[feature]])) +
    geom_violin(aes(fill = conservation_classification), alpha = 1, trim = TRUE) +
    scale_color_manual(values = c("tan1", "skyblue")) +
    scale_fill_manual(values = c("tan1", "skyblue")) +
    geom_signif(
      comparisons = list(c("cloud", "core")),
      map_signif_level = TRUE
    ) +
    theme_minimal() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(x = "", y = feature_titles[feature]) +
    ggtitle(feature_titles[feature])
})

# Save the arranged grid with updated titles
ggsave("Escherichia_coli_violin.png", 
       plot = do.call(grid.arrange, 
                      c(plots, 
                        ncol = 6, 
                        top = "Escherichia coli")), 
       width = 20, height = 14)


#### C BOTULINUM ####

df <- read_csv("/stor/scratch/Ochman/kristen/pangenome/all_bacterial_species/Clostridium_botulinum/rep_seq_properties/Clostridium_botulinum_seq_properties.csv", 
               col_names = FALSE)
colnames(df) <- c("gene_id", "feature", "value")

wide_df <- df %>% 
  pivot_wider(names_from = feature, values_from = value) %>% # 16756 rows
  dplyr::filter((conservation_percentage >= 0.95 & conservation_percentage <= 1) | conservation_percentage <= 0.05) %>% # 14043 rows
  mutate(conservation_classification = ifelse(conservation_percentage >= 0.95, "core", "cloud")) %>% 
  mutate(transmembrane = ifelse(transmembrane == 0, NA, transmembrane)) %>% 
  mutate(tm_coverage = ifelse(tm_coverage == 0, NA, tm_coverage)) %>% 
  mutate(across(all_of(features), 
                ~ ifelse(. >= quantile(., 0.05, na.rm = TRUE) & . <= quantile(., 0.95, na.rm = TRUE), ., NA)))

plots <- lapply(features, function(feature) {
  ggplot(wide_df, aes(x = conservation_classification, y = .data[[feature]])) +
    geom_violin(aes(fill = conservation_classification), alpha = 0.5, trim = TRUE) +
    scale_color_manual(values = c("tan1", "skyblue")) +
    scale_fill_manual(values = c("tan1", "skyblue")) +
    geom_signif(
      comparisons = list(c("cloud", "core")),
      map_signif_level = TRUE
    ) +
    theme_minimal() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(x = "", y = feature_titles[feature]) +
    ggtitle(feature_titles[feature])
})

# Save the arranged grid with updated titles
ggsave("Clostridium_botulinum_violin.png", 
       plot = do.call(grid.arrange, 
                      c(plots, 
                        ncol = 6, 
                        top = "Clostridium botulinum")), 
       width = 20, height = 14)
