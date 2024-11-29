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
        tm_value = parts[1]
        writer.writerow([seq_id, 'transmembrane', tm_value])
