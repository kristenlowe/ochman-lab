import os
import re
from Bio import SeqIO

# Define the base directory containing all bacterial species folders
base_dir = "/stor/scratch/Ochman/kristen/pangenome/all_bacterial_species"

# Loop through each bacterial species directory
for species_dir in os.listdir(base_dir):
    species_path = os.path.join(base_dir, species_dir)
    
    # Paths to input and output files
    clustering_dir = os.path.join(species_path, f"{species_dir}_proteins/clustering")
    cds_file = os.path.join(species_path, f"{species_dir}_CDSs", f"all_{species_dir}_cds.fna")
    protein_file = os.path.join(species_path, f"{species_dir}_proteins", f"all_{species_dir}_proteins.faa")
    output_cds_file = os.path.join(species_path, f"rep_{species_dir}_cds.fna")
    output_protein_file = os.path.join(species_path, f"rep_{species_dir}_proteins.faa")
    
    # Check if the necessary input files exist
    if not (os.path.isfile(cds_file) and os.path.isfile(protein_file)):
        continue
    
    # Read representative sequence IDs from clusters.tsv
    clusters_file = os.path.join(clustering_dir, "clusters.tsv")
    if not os.path.isfile(clusters_file):
        continue
    
    rep_seqs = set()
    with open(clusters_file, 'r') as f:
        for line in f:
            rep_seq = line.strip().split()[0]  # First column is the representative sequence
            rep_seqs.add(rep_seq)
    
    # Define a function to filter, deduplicate, and rename sequences
    def filter_deduplicate_sequences(input_file, output_file, rep_seqs):
        unique_sequences = set()  # Track unique sequence IDs to avoid duplicates
        
        with open(output_file, 'w') as output_handle:
            for record in SeqIO.parse(input_file, "fasta"):
                # Extract the protein ID from the record ID (e.g., WP_011407161.1)
                protein_id_match = re.search(r"WP_\d+\.\d+", record.id)
                if protein_id_match:
                    protein_id = protein_id_match.group()
                    # Check if the protein ID is a representative sequence and not already written
                    if protein_id in rep_seqs and protein_id not in unique_sequences:
                        # Add to unique sequences set to prevent duplicates
                        unique_sequences.add(protein_id)
                        # Rename the sequence ID and write it to the output file
                        record.id = protein_id
                        record.description = ""  # Clear the description
                        SeqIO.write(record, output_handle, "fasta")
    
    # Filter, deduplicate, and rename sequences for CDS and protein files
    filter_deduplicate_sequences(cds_file, output_cds_file, rep_seqs)
    filter_deduplicate_sequences(protein_file, output_protein_file, rep_seqs)

