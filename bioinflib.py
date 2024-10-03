#! /usr/bin/env python

from typing import List, Tuple, Union, Dict
from bioinflib_folder_1 import RNA_DNA_Module as RD
from bioinflib_folder_1 import FILTER_Module as FILTER


def run_dna_rna_tools(*args: Union[str, int]) -> Union[str, List[str]]:
    """
    Выполняет различные операции с ДНК и РНК последовательностями.

    Args:
        *args: Аргументы включают:
            - последовательности ДНК или РНК (строки)
            - последнее значение должно быть строкой, обозначающей действие
              ('complement', 'reverse_complement', 'transcribe', 'reverse',
              'what_is_that', 'MW', 'GC', 'Tm')

    Returns:
        str или List[str]: Результат обработки последовательностей.
        Если введена одна последовательность, возвращается строка,
        если несколько — список из строк.
    """
    action = args[-1]
    empty_list = []
    for seq in args[0:-1]:
        if action == 'complement':
            if RD.is_dna(seq):
                complement_seq = RD.complement_DNA(seq)
                empty_list.append(complement_seq)
            elif RD.is_rna(seq):
                complement_seq = RD.complement_RNA(seq)
                empty_list.append(complement_seq)
            else:
                empty_list.append(f'{seq} is not RNA or DNA')

        elif action == 'reverse_complement':
            if RD.is_dna(seq):
                rev_complement_seq = RD.reverse_complement_dna(seq)
                empty_list.append(rev_complement_seq)
            elif RD.is_rna(seq):
                rev_complement_seq = RD.reverse_complement_rna(seq)
                empty_list.append(rev_complement_seq)
            else:
                empty_list.append(f'{seq} is not RNA or DNA')

        elif action == 'transcribe':
            if RD.is_dna(seq) or RD.is_rna(seq):
                trans_seq = RD.transcribe(seq)
                empty_list.append(trans_seq)
            else:
                empty_list.append(f'{seq} is not RNA or DNA')

        elif action == 'reverse':
            if RD.is_dna(seq) or RD.is_rna(seq):
                rev_seq = RD.reverse(seq)
                empty_list.append(rev_seq)
            else:
                empty_list.append(f'{seq} is not RNA or DNA')

        elif action == 'what_is_that':
            if RD.is_dna(seq):
                empty_list.append(f'{seq} is DNA')
            elif RD.is_rna(seq):
                empty_list.append(f'{seq} is RNA')
            else:
                empty_list.append(f"I don't know what {seq} is")

        elif action == 'MW':
            if RD.is_dna(seq):
                MW_seq = RD.ssDNA_MW(seq)
                empty_list.append(f'{MW_seq} Da for this DNA')
            elif RD.is_rna(seq):
                MW_seq = RD.ssRNA_MW(seq)
                empty_list.append(f'{MW_seq} Da for this RNA')
            else:
                empty_list.append(f"I don't know what {seq} is")

        elif action == 'GC':
            if RD.is_dna(seq) or RD.is_rna(seq):
                GC_seq = RD.GC_content(seq)
                empty_list.append(f'GC content is {GC_seq} %')
            else:
                empty_list.append(f"I don't know what {seq} is")

        elif action == 'Tm':
            if RD.is_dna(seq):
                Tm_seq = RD.Tm_primer(seq)
                empty_list.append(f'Tm is {Tm_seq} degrees Celsius')
            else:
                empty_list.append(f"I can't do it for {seq}")

    if len(empty_list) > 1:
        return empty_list

    else:
        return empty_list[0]


def filter_fastq(seqs: Dict[str, Tuple[str, str]],
                 gc_bounds: Union[Tuple[float, float], float, int] = (0, 100),
                 length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
                 quality_threshold: int = 0) -> Dict[str, Tuple[str, str]]:
    """
    Отбирает последовательности в формате fastq по различным параметрам.

    Args:
        seqs: Словарь с fastq-ридами. Ключ - имя последовательности,
              значение - кортеж из двух строк: сама последовательность и ее качество.
        gc_bounds: Интервал GC состава (в процентах) для фильтрации.
                   Можно передать одно значение для указания верхней границы.
        length_bounds: Интервал длины последовательности для фильтрации.
                       Можно передать одно значение для указания верхней границы.
        quality_threshold: Пороговое значение среднего качества рида для фильтрации.
                           По-умолчанию 0 (шкала phred33).

    Returns:
        Dict[str, Tuple[str, str]]: Словарь с отфильтрованными последовательностями.
    """
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)

    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, length_bounds)

    filtered_seqs = {}

    for name, (sequence, quality) in seqs.items():

        seq_length = len(sequence)
        if length_bounds[0] <= seq_length <= length_bounds[1]:

            gc_percentage = FILTER.gc_content(sequence)
            if gc_bounds[0] <= gc_percentage <= gc_bounds[1]:

                avg_quality = FILTER.average_quality(quality)
                if avg_quality >= quality_threshold:

                    filtered_seqs[name] = (sequence, quality)

    return filtered_seqs
