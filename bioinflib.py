from typing import List, Tuple, Union, Dict
from bioinflib_folder_1 import RNA_DNA_Module as RD
from bioinflib_folder_1 import FILTER_Module as FILTER
from bioinflib_folder_1 import FASTQ_Module as FQ


def run_dna_rna_tools(*args: Union[str, int]) -> Union[str, List[str]]:
    """
    Performs various operations with DNA and RNA sequences.

    Args:
        *args: Arguments include:
            - DNA or RNA sequences (strings)
            - the last value should be a string indicating an action
              ('complement', 'reverse_complement', 'transcribe', 'reverse',
              'what_is_that', 'MW', 'GC', 'Tm')

    Returns:
        str or List[str]: The result of processing sequences.
        If one sequence is entered, the string is returned,
        if there are several, a list of strings.
    """
    action = args[-1]
    result = []
    for seq in args[0:-1]:
        if action == 'complement':
            if RD.is_dna(seq):
                complement_seq = RD.complement_DNA(seq)
                result.append(complement_seq)
            elif RD.is_rna(seq):
                complement_seq = RD.complement_RNA(seq)
                result.append(complement_seq)
            else:
                print(f'{seq} is not RNA or DNA')

        elif action == 'reverse_complement':
            if RD.is_dna(seq):
                rev_complement_seq = RD.reverse_complement_dna(seq)
                result.append(rev_complement_seq)
            elif RD.is_rna(seq):
                rev_complement_seq = RD.reverse_complement_rna(seq)
                result.append(rev_complement_seq)
            else:
                print(f'{seq} is not RNA or DNA')

        elif action == 'transcribe':
            if RD.is_dna(seq) or RD.is_rna(seq):
                trans_seq = RD.transcribe(seq)
                result.append(trans_seq)
            else:
                print(f'{seq} is not RNA or DNA')

        elif action == 'reverse':
            if RD.is_dna(seq) or RD.is_rna(seq):
                rev_seq = RD.reverse(seq)
                result.append(rev_seq)
            else:
                print(f'{seq} is not RNA or DNA')

        elif action == 'what_is_that':
            if RD.is_dna(seq):
                result.append(f'{seq} is DNA')
            elif RD.is_rna(seq):
                result.append(f'{seq} is RNA')
            else:
                result.append(f"I don't know what {seq} is")

        elif action == 'MW':
            if RD.is_dna(seq):
                MW_seq = RD.find_ssDNA_MW(seq)
                result.append(MW_seq)
            elif RD.is_rna(seq):
                MW_seq = RD.find_ssRNA_MW(seq)
                result.append(MW_seq)
            else:
                print(f"I don't know what {seq} is")

        elif action == 'GC':
            if RD.is_dna(seq) or RD.is_rna(seq):
                GC_seq = RD.find_GC_content(seq)
                result.append(f'{GC_seq}%')
            else:
                print(f"I don't know what {seq} is")

        elif action == 'Tm':
            if RD.is_dna(seq):
                Tm_seq = RD.find_Tm_primer(seq)
                result.append(f'{Tm_seq} degrees Celsius')
            else:
                print(f"I can't do it for {seq}")

    if len(result) > 1:
        return result

    else:
        return result[0]


def filter_fastq(input_fastq: str,
                 output_fastq: str,
                 gc_bounds: Union[Tuple[float, float], float, int] = (0, 100),
                 length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
                 quality_threshold: int = 0):
    """
    Selects sequences in fastq format according to various parameters.

    Args:
        input_fastq: The path to the input FASTQ file.
        output_fastq: The path to the output FASTQ file.
        gc_bounds: The GC composition interval (in percent) for filtering.
                   You can pass a single value to indicate the upper limit.
        length_bounds: The length interval of the sequence to filter.
                       You can pass a single value to indicate the upper limit.
        quality_threshold: The threshold value of the average read quality.
                           The default value is 0 (phred33 scale).

    Returns:
        Dict[str, Type[str, str]]: Dictionary with filtered sequences.
    """
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)

    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, length_bounds)

    for name, sequence, quality in FQ.read_fastq(input_fastq):
        seq_length = len(sequence)
        gc_percentage = FILTER.find_gc_content(sequence)
        avg_quality = FILTER.average_quality(quality)
        if ((length_bounds[0] <= seq_length <= length_bounds[1])
            and (gc_bounds[0] <= gc_percentage <= gc_bounds[1])
            and (avg_quality >= quality_threshold)):

            FQ.write_fastq(output_fastq, name, sequence, quality)
