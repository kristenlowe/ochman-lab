# generates a reversed sequence fasta file from a normal one
# (to make decoy sequences)

def reverse_sequences():
    # Prompt for input and output file paths
    input_file = input("Enter the input file path: ")
    output_file = input("Enter the output file path: ")

    # Read the input file
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Parse the sequences and reverse them
    sequences = {}
    current_seq = ''
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            current_seq = line
            sequences[current_seq] = ''
        else:
            sequences[current_seq] += line[::-1]

    # Write the reversed sequences to the output file
    with open(output_file, 'w') as f:
        for seq, rev_seq in sequences.items():
            f.write('{}\n{}\n'.format(seq, rev_seq))


# Call the function to reverse the sequences
reverse_sequences()
