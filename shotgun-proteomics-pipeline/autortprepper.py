import sys
import pandas as pd

def preprocess_peptide(peptide):
    # Remove modification information for C and replace modified M with 1
    peptide = peptide.replace('+57.021', '').replace('M+15.995', '1')
    return peptide

def main():
    # Check if the correct number of command line arguments are provided
    if len(sys.argv) != 5:
        print("Usage: python3 autortprepper.py input_msgf_file.tsv input_retention_file.txt test_file.tsv training_file.tsv")
        sys.exit(1)

    # Read input files
    input_msgf_file = sys.argv[1]
    input_retention_file = sys.argv[2]
    test_output_file = sys.argv[3]
    training_output_file = sys.argv[4]

    # Read the MSGF file
    msgf_df = pd.read_csv(input_msgf_file, sep='\t', comment='#')

    # Filter rows based on Q-value
    filtered_df = msgf_df[msgf_df.iloc[:, 15] < 0.01]

    # Ensure we don't exceed the available count
    test_count = min(1000, len(filtered_df))
    training_count = max(0, len(filtered_df) - test_count)

    # Extract peptides and preprocess
    test_peptides = filtered_df.iloc[:test_count, 9].apply(preprocess_peptide).tolist()

    # Read the retention time file
    retention_df = pd.read_csv(input_retention_file, sep='\t', header=None, names=['Peptide', 'RetentionTime'])

    # Filter rows based on test peptides
    test_retention_df = retention_df[retention_df['Peptide'].isin(test_peptides)]

    # Write test file
    test_retention_df.to_csv(test_output_file, sep='\t', columns=['Peptide', 'RetentionTime'], header=['x','y'], index=False)

    # Filter rows based on training peptides
    training_retention_df = retention_df[~retention_df['Peptide'].isin(test_peptides)].head(training_count)

    # Write training file
    training_retention_df.to_csv(training_output_file, sep='\t', columns=['Peptide', 'RetentionTime'], header=['x','y'], index=False)

if __name__ == "__main__":
    main()
