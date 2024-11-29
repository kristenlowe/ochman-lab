import csv
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    reader = infile.readlines()
    writer = csv.writer(outfile)

    for line in reader[1:]:  # Skip the header line
        parts = line.split()
        seq_id = parts[0]
        sp_value = parts[2]
        if sp_value == "Y":
            writer.writerow([seq_id, 'signal_peptide', 1])
        else:
            writer.writerow([seq_id, 'signal_peptide', sp_value])
