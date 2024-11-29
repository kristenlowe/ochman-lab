import csv
import pandas as pd
from tqdm import tqdm

# Step 1: Load the clusters.tsv file and create a reverse lookup dictionary
gene_to_family = {}
with open('/stor/scratch/Ochman/kristen/pangenome/seqprop_pipeline/all_419_CDS.NR.clusters.csv', 'r') as f:
    reader = csv.reader(f, delimiter=',')
    for row in reader:
        gene_family, gene = row
        gene_to_family[gene] = gene_family

# Step 2: Load the presence_absence_matrix.csv file and determine the number of strains each gene family is present in
presence_absence = pd.read_csv('/stor/scratch/Ochman/kristen/pangenome/seqprop_pipeline/all_419_CDS.NR.clusters.presenceabsence.csv', index_col=0)
gene_family_strain_count = {}

for gene_family in tqdm(presence_absence.index, desc="Processing presence/absence matrix"):
    num_strains = presence_absence.loc[gene_family].sum()
    gene_family_strain_count[gene_family] = num_strains

# Step 3: Load the ECOR_sequenceproperties.csv file and update it with the number of strains
updated_rows = []
with open('cleaned_all_419_seq_properties.csv', 'r') as f:
    reader = csv.reader(f)
    for row in tqdm(reader, desc="Updating sequence properties"):
        gene, property_type, value = row
        gene_family = gene_to_family.get(gene)
        if gene_family:
            num_strains = gene_family_strain_count[gene_family]
            updated_rows.append([gene, property_type, value, num_strains])

# Step 4: Write the updated information to a new CSV file without a header
with open('cleaned_all_419_seq_properties_w_freq.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(updated_rows)

print("Updated file has been written to 'cleaned_all_419_seq_properties_w_freq.csv'.")
