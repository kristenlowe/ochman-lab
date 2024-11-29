library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(gridExtra)
library(grid)
library(forcats)
library(ggsignif)
library(ggrepel)

setwd("/stor/scratch/Ochman/kristen/pangenome")

# import data
lifestyles <- read_csv("57_pangenome_lifestyles.csv")
species_pan <- read_csv("bac_57_difference_levels.csv")

life_df <- inner_join(x = species_pan, y = lifestyles, 
                      by = c("species" = "Species"))

#### HISTOGRAMS ####
create_overlapping_histogram <- function(data, property_filter, x_label, title, binwidth = NULL) {
  data %>% 
    pivot_longer(cols = c(core_value, cloud_value), 
                 names_to = "Group", 
                 values_to = "Value") %>% 
    mutate(Group = ifelse(Group == "core_value", "core", "cloud")) %>% 
    filter(property == property_filter) %>% 
    ggplot(aes(x = Value, fill = Group)) +
    geom_histogram(color = 'white', binwidth = binwidth) +
    scale_fill_manual(values = c("core" = "skyblue", "cloud" = "tan1")) +
    theme_bw() +
    labs(y = 'Number of Species',
         x = x_label,
         title = title)
}

# Length
create_overlapping_histogram(
  data = life_df,
  property_filter = 'length',
  x_label = 'Length (bp)',
  title = 'Distribution of Lengths',
  binwidth = 50
)

# GC
create_overlapping_histogram(
  data = life_df,
  property_filter = 'GC',
  x_label = 'GC Content (%)',
  title = 'Distribution of GC Content',
  binwidth = 5
)

# GC_3rd
create_overlapping_histogram(
  data = life_df,
  property_filter = 'GC_3rd',
  x_label = 'GC Content at 3rd Codon Position (%)',
  title = 'Distribution of GC Content at 3rd Codon Position',
  binwidth = 5
)

# Polar Amino Acid Content
create_overlapping_histogram(
  data = life_df,
  property_filter = 'polarAA',
  x_label = 'Polar Amino Acid Content (%)',
  title = 'Distribution of Polar Amino Acid Content',
  binwidth = 2
)

# Hydrophobic Amino Acid Content
create_overlapping_histogram(
  data = life_df,
  property_filter = 'hydrophobicAA',
  x_label = 'Hydrophobic Amino Acid Content (%)',
  title = 'Distribution of Hydrophobic Amino Acid Content',
  binwidth = 2
)

# Metabolic Cost
create_overlapping_histogram(
  data = life_df,
  property_filter = 'metabol',
  x_label = 'Metabolic Cost (ATP/mol AA)',
  title = 'Distribution of Metabolic Cost',
  binwidth = 0.3
)

# Codon Adaptation Index
create_overlapping_histogram(
  data = life_df,
  property_filter = 'cai',
  x_label = 'Codon Adaptation Index',
  title = 'Distribution of Codon Adaptation Index',
  binwidth = 0.03
)

# Instability Index
create_overlapping_histogram(
  data = life_df,
  property_filter = 'instability',
  x_label = 'Instability Index',
  title = 'Distribution of Instability Index',
  binwidth = 1
)

# Intrinsic Structural Disorder
create_overlapping_histogram(
  data = life_df,
  property_filter = 'disorder',
  x_label = 'Intrinsic Structural Disorder',
  title = 'Distribution of Intrinsic Structural Disorder',
  binwidth = 0.02
)

# Transmembrane Domain Count
create_overlapping_histogram(
  data = life_df,
  property_filter = 'transmembrane',
  x_label = 'Transmembrane Domain Count',
  title = 'Distribution of Transmembrane Domain Count',
  binwidth = 1
)

# Transmembrane Coverage
create_overlapping_histogram(
  data = life_df,
  property_filter = 'tm_coverage',
  x_label = 'Transmembrane Coverage',
  title = 'Distribution of Transmembrane Coverage',
  binwidth = 0.05
)

# Amino Acid Composition Bias
create_overlapping_histogram(
  data = life_df,
  property_filter = 'aa_comp_bias',
  x_label = 'Amino Acid Composition Bias',
  title = 'Distribution of Amino Acid Composition Bias',
  binwidth = 1
)


#### SCATTERPLOTS ####

# Create a new data frame
life_forscatter <- life_df %>%
  group_by(species) %>%
  mutate(core_gc_content = core_value[property == "GC"]) %>%
  ungroup()

# Generate random numbers between 5 and 290 for rows where p-value is 0
random_exponents <- runif(
  sum(life_forscatter$pvalue == 0), 
  min = 5, 
  max = 290
)

# Replace p-values equal to 0 with 1e-(random number)
life_forscatter$pvalue[life_forscatter$pvalue == 0] <- 10^(-random_exponents)

create_scatterplot <- function(data, 
                               property_filter, 
                               x_label, 
                               title, 
                               legend_label) {
  data %>% 
    filter(property == property_filter) %>% 
    ggplot(aes(x = core_minus_cloud_difference, 
               y = -log10(pvalue),
               color = core_gc_content)) +
    geom_jitter(alpha = 0.7) +
    theme_bw() +
    geom_hline(yintercept = -log10(0.05), 
               linetype = "dashed", 
               color = "red") +
    geom_vline(xintercept = 0, 
               linetype = "dashed", 
               color = "black") +
    labs(y = "-log10(p-value) (Significance)",
         x = x_label,
         title = title,
         color = "GC Content (%)")
}


# Length
create_scatterplot(
  data = life_forscatter,
  property_filter = 'length',
  x_label = "Difference in Gene Lengths (Core - Cloud) (bp)",
  title = 'Length Differences and Significance',
  legend_label = 'Core Gene\nLength (bp)'
) +
  geom_text_repel(data = life_forscatter %>% 
                    filter(species == 'Mycobacterium_tuberculosis') %>% 
                    filter(property == 'length'), 
                  aes(label = species), 
                  size = 3)

# GC
create_scatterplot(
  data = life_forscatter,
  property_filter = 'GC',
  x_label = "Difference in GC Content (Core - Cloud) (%)",
  title = 'GC Content Differences and Significance',
  legend_label = 'Core GC\nContent (%)'
) +
  geom_text_repel(data = life_forscatter %>% 
                    filter(species == 'Mycobacterium_tuberculosis' | 
                             species == 'Aeromonas_hydrophila' | 
                             species == 'Chlamydia_trachomatis') %>% 
                    filter(property == 'GC'), 
                  aes(label = species), 
                  size = 3)

# GC_3rd
create_scatterplot(
  data = life_forscatter,
  property_filter = 'GC_3rd',
  x_label = "Difference in GC Content at 3rd Codon Position (Core - Cloud) (%)",
  title = 'GC Content at 3rd Codon Position Differences and Significance',
  legend_label = 'Core GC\nContent at\n3rd Codon\nPosition (%)'
) +
  geom_text_repel(data = life_forscatter %>% 
                    filter(species == 'Chlamydia_trachomatis' |
                             species == 'Streptococcus_agalactiae' |
                             species == 'Clostridioides_difficile' |
                             species == 'Aeromonas_hydrophila' |
                             species == 'Serratia_marcescens' |
                             species == 'Klebsiella_pneumoniae') %>% 
                    filter(property == 'GC_3rd'), 
                  aes(label = species), 
                  size = 3)

# polarAA
create_scatterplot(
  data = life_forscatter,
  property_filter = 'polarAA',
  x_label = "Difference in Polar Amino Acid Content (Core - Cloud) (%)",
  title = 'Polar Amino Acid Content Differences and Significance',
  legend_label = 'Core Polar\nAmino Acid\nContent (%)'
) +
  geom_text_repel(data = life_forscatter %>% 
                    filter(species == 'Mycobacterium_tuberculosis') %>% 
                    filter(property == 'polarAA'), 
                  aes(label = species), 
                  size = 3)

# hydrophobicAA
create_scatterplot(
  data = life_forscatter,
  property_filter = 'hydrophobicAA',
  x_label = "Difference in Hydrophobic Amino Acid Content (Core - Cloud) (%)",
  title = 'Hydrophobic Amino Acid Content Differences and Significance',
  legend_label = 'Core Hydrophobic\nAmino Acid\nContent (%)'
) +
  geom_text_repel(data = life_forscatter %>% 
                    filter(species == 'Mycobacterium_tuberculosis') %>% 
                    filter(property == 'hydrophobicAA'), 
                  aes(label = species), 
                  size = 3)

# metabol
create_scatterplot(
  data = life_forscatter,
  property_filter = 'metabol',
  x_label = "Difference in Metabolic Cost (Core - Cloud) (ATP/mol AA)",
  title = 'Metabolic Cost Differences and Significance',
  legend_label = 'Core Metabolic\nCost (ATP/mol AA)'
) +
  geom_text_repel(data = life_forscatter %>% 
                    filter(species == 'Mycobacterium_tuberculosis') %>% 
                    filter(property == 'metabol'), 
                  aes(label = species), 
                  size = 3)

# cai
create_scatterplot(
  data = life_forscatter,
  property_filter = 'cai',
  x_label = "Difference in Codon Adaptation Index (Core - Cloud)",
  title = 'Codon Adaptation Index Differences and Significance',
  legend_label = 'Core Codon\nAdaptation\nIndex'
) +
  geom_text_repel(data = life_forscatter %>% 
                    filter(species == 'Mycobacterium_tuberculosis' |
                             species == 'Clostridium_botulinum' |
                             species == 'Clostridioides_difficile') %>% 
                    filter(property == 'cai'), 
                  aes(label = species), 
                  size = 3)

# instability
create_scatterplot(
  data = life_forscatter,
  property_filter = 'instability',
  x_label = "Difference in Instability Index (Core - Cloud)",
  title = 'Instability Index Differences and Significance',
  legend_label = 'Core\nInstability\nIndex'
) +
  geom_text_repel(data = life_forscatter %>% 
                    filter(species == 'Mycobacterium_tuberculosis') %>% 
                    filter(property == 'instability'), 
                  aes(label = species), 
                  size = 3)

# disorder
create_scatterplot(
  data = life_forscatter,
  property_filter = 'disorder',
  x_label = "Difference in Instrinsic Structural Disorder (Core - Cloud)",
  title = 'Instrinsic Structural Disorder Differences and Significance',
  legend_label = 'Core Intrinsic\nStructural Disorder'
) +
  geom_text_repel(data = life_forscatter %>% 
                    filter(species == 'Mycobacterium_tuberculosis' |
                             species == 'Clostridium_botulinum') %>% 
                    filter(property == 'disorder'), 
                  aes(label = species), 
                  size = 3)

# transmembrane
create_scatterplot(
  data = life_forscatter,
  property_filter = 'transmembrane',
  x_label = "Difference in Transmembrane Domain Count (Core - Cloud)",
  title = 'Transmembrane Domain Count Differences and Significance',
  legend_label = 'Core Transmembrane\nDomain Count'
) +
  geom_text_repel(data = life_forscatter %>% 
                    filter(species == 'Staphylococcus_aureus') %>% 
                    filter(property == 'transmembrane'), 
                  aes(label = species), 
                  size = 3)

# tm_coverage
create_scatterplot(
  data = life_forscatter,
  property_filter = 'tm_coverage',
  x_label = "Difference in Transmembrane Coverage (Core - Cloud)",
  title = 'Transmembrane Coverage Differences and Significance',
  legend_label = 'Core Transmembrane\nCoverage'
) +
  geom_text_repel(data = life_forscatter %>% 
                    filter(species == 'Staphylococcus_aureus') %>% 
                    filter(property == 'tm_coverage'), 
                  aes(label = species), 
                  size = 3)

# aa_comp_bias
create_scatterplot(
  data = life_forscatter,
  property_filter = 'aa_comp_bias',
  x_label = "Difference in Amino Acid Composition Bias (Core - Cloud)",
  title = 'Amino Acid Composition Bias Differences and Significance',
  legend_label = 'Core Amino Acid\nComposition Bias'
) +
  geom_text_repel(data = life_forscatter %>% 
                    filter(species == 'Corynebacterium_pseudotuberculosis' |
                             species == 'Clostridium_botulinum') %>% 
                    filter(property == 'aa_comp_bias'), 
                  aes(label = species), 
                  size = 3)








