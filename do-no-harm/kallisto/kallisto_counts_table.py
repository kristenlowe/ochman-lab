#!/usr/bin/env python3

import csv
import os
from collections import defaultdict

# Define the array of kmers
kmers = [7, 11, 15, 19, 23, 27, 31]

# Process the abundance files and create the counts table
for kmer in kmers:
    output_filename = f'run{kmer}_counts_table.csv'
    sample_counts = defaultdict(lambda: defaultdict(int))  # Nested dictionary for counts
    samples = []

    # Search for directories and process files
    for dir_entry in os.scandir():
        if dir_entry.is_dir() and f'_run{kmer}' in dir_entry.name:
            sample_name = dir_entry.name.split('_run')[0]
            samples.append(sample_name)
            file_path = os.path.join(dir_entry.path, 'abundance.tsv')

            with open(file_path, 'r') as tsvfile:
                reader = csv.DictReader(tsvfile, delimiter='\t')
                for row in reader:
                    target_id = row['target_id']
                    est_counts = float(row['est_counts'])
                    sample_counts[target_id][sample_name] = est_counts

    # Write to the output CSV file
    with open(output_filename, 'w', newline='') as csvfile:
        fieldnames = ['target_id'] + samples
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for target_id in sample_counts:
            row = {'target_id': target_id}
            for sample in samples:
                row[sample] = sample_counts[target_id].get(sample, 0)
            writer.writerow(row)
