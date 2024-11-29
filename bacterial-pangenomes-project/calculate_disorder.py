import csv
from tqdm import tqdm
import sys

def parse_fasta_file(fasta_path):
    sequences = {}
    with open(fasta_path, 'r') as file:
        current_header = None
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_header = line.split()[0][1:]
                sequences[current_header] = []
            else:
                sequences[current_header].append(line)
    for key in sequences:
        sequences[key] = ''.join(sequences[key])
    return sequences

def parse_iupred_output(iupred_path):
    with open(iupred_path, 'r') as file:
        data = file.read().splitlines()

    scores = []
    for line in tqdm(data, desc="Processing IUPred output"):
        if line and line[0].isdigit():
            parts = line.split()
            if len(parts) == 3:
                score = float(parts[2])
                scores.append(score)
    return scores

def calculate_average_disorder(sequences, scores):
    average_disorder = {}
    score_index = 0
    for protein, sequence in tqdm(sequences.items(), desc="Calculating averages"):
        protein_scores = []
        for residue in sequence:
            if score_index < len(scores):
                protein_scores.append(scores[score_index])
                score_index += 1
            else:
                print(f"Warning: Ran out of scores before matching all residues of {protein}")
                break
        avg_score = sum(protein_scores) / len(protein_scores) if protein_scores else 0
        average_disorder[protein] = avg_score
    return average_disorder

def write_to_csv(output_path, average_disorder_scores):
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for protein, avg_score in tqdm(average_disorder_scores.items(), desc="Writing to CSV"):
            writer.writerow([protein, 'disorder', avg_score])

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python calculate_disorder.py <fasta_file_path> <iupred_output_path> <output_csv_path>")
        sys.exit(1)

    fasta_file_path = sys.argv[1]
    iupred_output_path = sys.argv[2]
    output_csv_path = sys.argv[3]

    sequences = parse_fasta_file(fasta_file_path)
    print(f"Parsed {len(sequences)} sequences from FASTA file.")

    iupred_scores = parse_iupred_output(iupred_output_path)
    print(f"Parsed {len(iupred_scores)} scores from IUPred output.")

    average_disorder_scores = calculate_average_disorder(sequences, iupred_scores)
    print(f"Calculated average disorder scores for {len(average_disorder_scores)} proteins.")

    write_to_csv(output_csv_path, average_disorder_scores)
    print(f"Results written to {output_csv_path}")

