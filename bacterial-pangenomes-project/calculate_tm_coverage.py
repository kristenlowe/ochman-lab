import csv
import sys

def parse_phobius_output(phobius_file):
    tm_data = {}
    with open(phobius_file, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            parts = line.split()
            sequence_id = parts[0]
            tm_count = int(parts[1])
            if tm_count == 0:
                tm_data[sequence_id] = 0
                continue

            prediction = parts[3]
            tm_domains = []
            segments = prediction.split('i')
            for segment in segments:
                if 'o' in segment:
                    domains = segment.split('o')
                    for domain in domains:
                        if '-' in domain and '/' not in domain:
                            start, end = map(int, domain.split('-'))
                            tm_domains.append((start, end))
            tm_data[sequence_id] = tm_domains
    return tm_data

def parse_lengths_file(lengths_file):
    lengths_data = {}
    with open(lengths_file, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            sequence_id = row[0]
            sequence_length = int(row[2])
            lengths_data[sequence_id] = sequence_length
    return lengths_data

def calculate_coverage(tm_data, lengths_data):
    coverage_data = []
    for sequence_id, sequence_length in lengths_data.items():
        if sequence_id in tm_data:
            if tm_data[sequence_id] == 0:
                coverage_data.append((sequence_id, 'tm_coverage', 0.0))
                continue
            total_tm_length = sum(end - start + 1 for start, end in tm_data[sequence_id])
            codon_length = sequence_length // 3
            coverage = total_tm_length / codon_length
            coverage_data.append((sequence_id, 'tm_coverage', coverage))
        else:
            coverage_data.append((sequence_id, 'tm_coverage', 0.0))
    return coverage_data

def write_coverage_to_csv(coverage_data, output_file):
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(coverage_data)

if __name__ == "__main__":
    phobius_file = sys.argv[1]
    lengths_file = sys.argv[2]
    output_file = sys.argv[3]
    
    tm_data = parse_phobius_output(phobius_file)
    lengths_data = parse_lengths_file(lengths_file)
    coverage_data = calculate_coverage(tm_data, lengths_data)
    write_coverage_to_csv(coverage_data, output_file)

    print(f"Coverage data written to {output_file}")
