import pandas as pd
from glob import glob

# Load cluster data
clusters = pd.read_csv("clusters.tsv", sep="\t", header=None, names=["protein1", "protein2"])

# Extract unique protein identifiers as gene families and unique strain identifiers
clusters['family'] = clusters['protein1']
clusters['strain1'] = clusters['protein1'].apply(lambda x: x.split('_')[1])
clusters['strain2'] = clusters['protein2'].apply(lambda x: x.split('_')[1])

# Create a set of all strains involved and convert it to a sorted list
all_strains = sorted(set(clusters['strain1']).union(set(clusters['strain2'])))

# Convert set of gene families to a sorted list
gene_families = sorted(set(clusters['family']))

# Initialize the matrix with zeros
presence_absence_matrix = pd.DataFrame(0, index=gene_families, columns=all_strains)

# Fill the matrix by marking presence
for index, row in clusters.iterrows():
    presence_absence_matrix.at[row['family'], row['strain1']] = 1
    presence_absence_matrix.at[row['family'], row['strain2']] = 1

# Save the matrix to a CSV file
presence_absence_matrix.to_csv("presence_absence_matrix.csv")
