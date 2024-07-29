import sys
import re

def replace_modified_m(sequence, modification_info):
    modified_positions = [int(match.group(1)) for match in re.finditer(r'(\d+)-\[Oxidation', modification_info)]
    sequence_list = list(sequence)
    for position in modified_positions:
        sequence_list[position - 1] = '1'
    return ''.join(sequence_list)

def filter_retention_times(peptide_data):
    filtered_data = {}
    for sequence, rt_values in peptide_data.items():
        if len(rt_values) == 1:
            filtered_data[sequence] = rt_values[0]
        else:
            rt_values.sort()
            min_rt = rt_values[0]
            max_rt = rt_values[-1]
            if max_rt - min_rt <= 180:
                avg_rt = sum(rt_values) / len(rt_values)
                filtered_data[sequence] = avg_rt
    return filtered_data

def retention_data(input_filename, output_filename):
    with open(input_filename, 'r') as input_file:
        lines = input_file.readlines()
        i = 0
        peptide_data = {}

        while i < len(lines):
            if lines[i].startswith("SEQ="):
                sequence = lines[i].split('=')[1].strip()
                i += 1
                while i < len(lines) and not lines[i].startswith("RTINSECONDS="):
                    if lines[i].startswith("USER03="):
                        modification_info = lines[i].split('=')[1].strip()
                    i += 1
                if i < len(lines) and lines[i].startswith("RTINSECONDS="):
                    rt_seconds = float(lines[i].split('=')[1].strip())
                    sequence = replace_modified_m(sequence, modification_info)
                    
                    if sequence not in peptide_data:
                        peptide_data[sequence] = [rt_seconds]
                    else:
                        peptide_data[sequence].append(rt_seconds)
            i += 1

    filtered_peptide_data = filter_retention_times(peptide_data)

    with open(output_filename, 'w') as output_file:
        for sequence, rt_value in filtered_peptide_data.items():
            output_file.write(f"{sequence}\t{rt_value}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 retention_data_extractor.py input_file.mgf output_file.txt")
    else:
        input_filename = sys.argv[1]
        output_filename = sys.argv[2]
        retention_data(input_filename, output_filename)
