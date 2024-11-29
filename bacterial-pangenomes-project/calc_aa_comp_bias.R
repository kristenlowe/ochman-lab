#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
protein_fasta_file <- args[1]
output_csv_file <- args[2]

# Load required library
library(Biostrings)

# Read the protein sequences
rep_proteins <- readAAStringSet(protein_fasta_file)

# Function to calculate reference composition
calculate_reference_composition <- function(proteins) {
  aa_counts <- alphabetFrequency(proteins, as.prob = TRUE)
  avg_composition <- colMeans(aa_counts)
  return(avg_composition)
}

# Function to calculate standard deviation of composition
calculate_sd_composition <- function(proteins, reference_composition) {
  compositions <- sapply(proteins, function(protein) {
    alphabetFrequency(protein, as.prob = TRUE)
  })
  sd_composition <- apply(compositions, 1, sd)
  return(sd_composition)
}

# Function to calculate composition bias
calculate_composition_bias <- function(protein, reference_composition, sd_composition) {
  protein_composition <- alphabetFrequency(protein, as.prob = TRUE)
  bias <- sum(abs(protein_composition - reference_composition) / sd_composition, na.rm = TRUE)
  return(bias)
}

# Calculate reference composition and standard deviation
reference_composition <- calculate_reference_composition(rep_proteins)
sd_composition <- calculate_sd_composition(rep_proteins, reference_composition)

# Calculate composition bias for each protein
composition_biases <- sapply(rep_proteins, calculate_composition_bias, reference_composition, sd_composition)

# Write the results to a CSV file
output_data <- data.frame(
  protein_id = names(rep_proteins),
  property = "aa_comp_bias",
  value = composition_biases
)

write.table(output_data, file = output_csv_file, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

