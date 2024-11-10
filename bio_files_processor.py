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


def parse_blast_output(input_file: str, output_file: str = None) -> None:
    """
    Parses the BLAST output, extracting a description of the first match for each request.
    
    Arguments:
    - input_file (str): Path to the input file with BLAST output.
    - output_file (str, optional): The path to the output file for recording match descriptions.
    """
    if output_file is None:
        output_file = f"{input_file.rsplit('.', 1)[0]}_result.txt"

    matches = []
    with open(input_file, 'r') as infile:
        for line in infile:
            line = line.strip()
            if 'Sequences producing significant alignments:' in line:
                infile.readline()
                infile.readline()
                next_line = infile.readline().strip()
                if next_line:
                    description = next_line.split('[', 1)[0]
                    matches.append(description)
    matches = sorted(matches)
    with open(output_file, 'w') as outfile:
        for hit in matches:
            outfile.write(f"{hit}\n")
