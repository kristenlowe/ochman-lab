# Load required libraries
library(tidyverse)

# Define the base directory containing the species folders and the output directory for PNGs
base_dir <- "/stor/scratch/Ochman/kristen/pangenome/all_bacterial_species"
output_dir <- "/stor/scratch/Ochman/kristen/pangenome"

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the features to plot
features <- c("length", "GC", "GC_3rd", "polarAA", "hydrophobicAA", "metabol", 
              "cai", "instability", "disorder", "aa_comp_bias", "transmembrane", "tm_coverage")

# Create an empty data frame to store the results
all_species_summary <- data.frame()

# Loop through each species directory and process the seq_properties.csv file
species_dirs <- list.dirs(base_dir, recursive = FALSE)

for (species_dir in species_dirs) {
  # Extract genus_species name from the directory path
  genus_species <- basename(species_dir)
  
  # Define the path to the seq_properties.csv file
  seq_prop_file <- file.path(species_dir, "rep_seq_properties", paste0(genus_species, "_seq_properties.csv"))
  
  # Check if the seq_properties.csv file exists
  if (file.exists(seq_prop_file)) {
    # Read the seq_properties.csv file
    df <- read_csv(seq_prop_file, col_names = FALSE)
    colnames(df) <- c("gene_id", "feature", "value")
    
    # Process data: pivot wider, filter, and classify conservation
    wide_df <- df %>% 
      pivot_wider(names_from = feature, values_from = value) %>%
      dplyr::filter((conservation_percentage >= 0.95 & conservation_percentage <= 1) | conservation_percentage <= 0.05) %>%
      mutate(conservation_classification = ifelse(conservation_percentage >= 0.95, "core", "cloud")) %>%
      mutate(transmembrane = ifelse(transmembrane == 0, NA, transmembrane)) %>% 
      mutate(tm_coverage = ifelse(tm_coverage == 0, NA, tm_coverage))
    
    # Compute the median values for cloud
    cloud_temp <- wide_df %>% 
      filter(conservation_classification == "cloud") %>% 
      dplyr::select(length, GC, GC_3rd, polarAA, hydrophobicAA, metabol, cai, 
                    instability, disorder, aa_comp_bias, transmembrane, tm_coverage) %>% 
      summarise(across(where(is.numeric), median, na.rm = TRUE)) %>% 
      mutate(conservation_classification = "cloud", species = genus_species)
    
    # Compute the median values for core
    core_temp <- wide_df %>% 
      filter(conservation_classification == "core") %>% 
      dplyr::select(length, GC, GC_3rd, polarAA, hydrophobicAA, metabol, cai, 
                    instability, disorder, aa_comp_bias, transmembrane, tm_coverage) %>% 
      summarise(across(where(is.numeric), median, na.rm = TRUE)) %>% 
      mutate(conservation_classification = "core", species = genus_species)
    
    # Append the results to the overall data frame
    all_species_summary <- bind_rows(all_species_summary, cloud_temp, core_temp)
    
  } else {
    message(paste("No seq_properties.csv found for", genus_species))
  }
}

# reorder columns
all_species_summary <- all_species_summary %>%
  select(species, conservation_classification, everything())

all_species_summary_long <- all_species_summary %>%
  pivot_longer(cols = length:tm_coverage, names_to = "property", values_to = "value") %>%
  pivot_wider(names_from = conservation_classification, values_from = value, names_prefix = "", values_fn = list(value = mean)) %>%
  rename(cloud_value = cloud, core_value = core) %>%
  select(species, property, cloud_value, core_value)


# Create empty columns for p-value and W test statistic in all_species_summary_long
all_species_summary_long <- all_species_summary_long %>%
  mutate(pvalue = NA_real_, W_test_statistic = NA_real_)

# Loop through each species and property to perform the Wilcoxon test
for (i in 1:nrow(all_species_summary_long)) {
  # Extract the species and property for the current row
  species <- all_species_summary_long$species[i]
  property <- all_species_summary_long$property[i]
  
  # Filter the original data (wide_df) for the current species
  species_dir <- file.path(base_dir, species)
  seq_prop_file <- file.path(species_dir, "rep_seq_properties", paste0(species, "_seq_properties.csv"))
  
  if (file.exists(seq_prop_file)) {
    # Read the seq_properties.csv file
    df <- read_csv(seq_prop_file, col_names = FALSE)
    colnames(df) <- c("gene_id", "feature", "value")
    
    # Process data: pivot wider, filter, and classify conservation
    wide_df <- df %>% 
      pivot_wider(names_from = feature, values_from = value) %>%
      dplyr::filter((conservation_percentage >= 0.95 & conservation_percentage <= 1) | conservation_percentage <= 0.05) %>%
      mutate(conservation_classification = ifelse(conservation_percentage >= 0.95, "core", "cloud")) %>%
      mutate(transmembrane = ifelse(transmembrane == 0, NA, transmembrane)) %>% 
      mutate(tm_coverage = ifelse(tm_coverage == 0, NA, tm_coverage))
    
    # Filter data for the core and cloud groups for the current property
    core_values <- wide_df %>%
      filter(conservation_classification == "core") %>%
      pull(!!sym(property))
    
    cloud_values <- wide_df %>%
      filter(conservation_classification == "cloud") %>%
      pull(!!sym(property))
    
    # Perform Wilcoxon test if both groups have values
    if (length(core_values) > 0 & length(cloud_values) > 0) {
      test_result <- wilcox.test(core_values, cloud_values)
      all_species_summary_long$pvalue[i] <- test_result$p.value
      all_species_summary_long$W_test_statistic[i] <- test_result$statistic
    }
  } else {
    message(paste("No seq_properties.csv found for", species))
  }
}

difference_df <- all_species_summary_long %>% 
  mutate(core_minus_cloud_difference = core_value - cloud_value)

write.csv(difference_df, file = "/stor/scratch/Ochman/kristen/pangenome/bac_57_difference_levels.csv", row.names = FALSE)


