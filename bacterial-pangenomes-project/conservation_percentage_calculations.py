import os
import re
import csv
from collections import defaultdict

# Define the base directory containing all bacterial species folders
base_dir = "/stor/scratch/Ochman/kristen/pangenome/all_bacterial_species"

# Loop through each bacterial species directory
for species_dir in os.listdir(base_dir):
    species_path = os.path.join(base_dir, species_dir)
    proteins_dir = os.path.join(species_path, f"{species_dir}_proteins")
    clustering_dir = os.path.join(proteins_dir, "clustering")
    clusters_file = os.path.join(clustering_dir, "clusters.tsv")
    
    # Check if the clusters file exists
    if not os.path.isfile(clusters_file):
        continue
    
    # Determine the number of strains by finding the highest strain number in the protein filenames
    num_strains = 0
    for filename in os.listdir(proteins_dir):
        match = re.search(rf"{species_dir}_(\d+)_protein\.faa", filename)
        if match:
            strain_num = int(match.group(1))
            num_strains = max(num_strains, strain_num)
    
    # Create the output directory and CSV file path
    rep_seq_properties_dir = os.path.join(species_path, "rep_seq_properties")
    os.makedirs(rep_seq_properties_dir, exist_ok=True)
    output_csv = os.path.join(rep_seq_properties_dir, f"{species_dir}_conservation.csv")
    
    # Initialize a dictionary to count sequences per representative
    cluster_counts = defaultdict(int)
    
    # Read each line in the clusters file and count occurrences of each representative
    with open(clusters_file, 'r') as f:
        for line in f:
            rep_seq, _ = line.strip().split()
            cluster_counts[rep_seq] += 1
    
    # Write the results to the CSV file
    with open(output_csv, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        for rep_seq, count in cluster_counts.items():
            conservation_percentage = count / num_strains
            csv_writer.writerow([rep_seq, "conservation_percentage", conservation_percentage])

