import csv
import sys
from Bio import SeqIO
from Bio.SeqUtils import ProtParamData
from tqdm import tqdm  # For the progress bar

def instability_index(seq):
    """Calculate the instability index according to Guruprasad et al 1990."""
    index = ProtParamData.DIWV  # Dipeptide instability weights
    score = 0.0
    length = len(seq)
    
    # Loop through each dipeptide and skip non-standard amino acids
    for i in range(length - 1):
        this, next = seq[i:i + 2]
        
        # Check if both amino acids are in the standard 20 amino acids
        if this in index and next in index[this]:
            dipeptide_value = index[this][next]
            score += dipeptide_value
        else:
            # Skip non-standard amino acids
            continue

    return (10.0 / length) * score

def calculate_instability_for_fasta(fasta_file, output_file):
    # Parse the FASTA file and count the total number of sequences for tqdm
    records = list(SeqIO.parse(fasta_file, "fasta"))
    
    # Open the output file for writing
    with open(output_file, mode='w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        # Write header
        #csvwriter.writerow(["gene", "feature", "value"])
        
        # Use tqdm to add a progress bar
        for record in tqdm(records, desc="Calculating Instability Index"):
            seq = str(record.seq)
            instability = instability_index(seq)
            csvwriter.writerow([record.id, "instability", instability])

# Specify the path to your FASTA file and the output CSV file
fasta_file = sys.argv[1]
output_file = sys.argv[2]

# Run the function to calculate instability and save to CSV
calculate_instability_for_fasta(fasta_file, output_file)

