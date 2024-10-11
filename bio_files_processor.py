def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None) -> None:
    """
    Converts a FASTA file with multiline sequences to a format where each sequence is written in one line.

    Arguments:
    input_data (str): The path to the input FASTA file with multiline sequences.
    output_fasta (str, optional): The path to the output FASTA file with sequences written in one line.
                                     If not specified, a file named "<input filename>_oneline.fasta".
    """
    if output_fasta is None:
        output_fasta = f"{input_fasta.rsplit('.', 1)[0]}_oneline.fasta"
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        header = ''
        sequence = ''
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    outfile.write(header + '\n')
                    outfile.write(sequence + '\n')
                header = line
                sequence = ''
            else:
                sequence += line
        if sequence:
            outfile.write(header + '\n')
            outfile.write(sequence + '\n')

# convert_multiline_fasta_to_oneline('example_multiline_fasta.fasta')
