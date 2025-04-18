import argparse
import logging
from typing import Union, Tuple
from Bio import SeqIO
from Bio import SeqUtils

def setup_logging(log_file: str = 'fastq_filter.log'):
    """Настройка логирования"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

def parse_gc_bounds(gc_arg: str) -> Union[Tuple[float, float], float]:
    """Парсинг аргумента gc_bounds"""
    if ',' in gc_arg:
        lower, upper = map(float, gc_arg.split(','))
        return (lower, upper)
    return float(gc_arg)

def parse_length_bounds(len_arg: str) -> Union[Tuple[int, int], int]:
    """Парсинг аргумента length_bounds"""
    if ',' in len_arg:
        lower, upper = map(int, len_arg.split(','))
        return (lower, upper)
    return int(len_arg)

def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: Union[float, int] = 0,
) -> int:
    """
    Function filter_fastq
    
    Args:
        input_fastq: The path to the input FASTQ file.
        output_fastq: The path to the output FASTQ file.
        gc_bounds: The GC composition interval (in percent) for filtering.
        length_bounds: The length interval of the sequence to filter.
        quality_threshold: The threshold value of the average read quality.
    
    Returns:
        int: Number of filtered records
    """
    try:
        with open(input_fastq, "rt") as in_handle:
            records = list(SeqIO.parse(in_handle, "fastq"))
            
            if not isinstance(gc_bounds, tuple):
                gc_bounds = (0, gc_bounds)
            if not isinstance(length_bounds, tuple):
                length_bounds = (0, length_bounds)

            logging.info(f"Starting filtering of {input_fastq}")
            logging.info(f"Parameters: GC={gc_bounds}, Length={length_bounds}, Quality={quality_threshold}")

            filtered = []
            for record in records:
                try:
                    qual = record.letter_annotations["phred_quality"]
                    mean_quality = sum(qual) / len(qual)
                    gc_content = SeqUtils.gc_fraction(record.seq) * 100
                    length = len(record.seq)

                    if (mean_quality >= quality_threshold and
                        gc_bounds[0] <= gc_content <= gc_bounds[1] and
                        length_bounds[0] <= length <= length_bounds[1]):
                        filtered.append(record)
                except Exception as e:
                    logging.warning(f"Skipping record {record.id}: {str(e)}")
                    continue

            with open(output_fastq, "wt") as out_handle:
                count = SeqIO.write(filtered, out_handle, "fastq")
            
            logging.info(f"Filtered {count} records to {output_fastq}")
            return count
            
    except FileNotFoundError as e:
        error_msg = f"Input file not found: {input_fastq}"
        logging.error(error_msg)
        raise FileNotFoundError(error_msg) from e
    except Exception as e:
        error_msg = f"Processing error: {str(e)}"
        logging.error(error_msg)
        raise RuntimeError(error_msg) from e

def main():
    parser = argparse.ArgumentParser(description='Filter FASTQ files by GC content, length and quality.')
    parser.add_argument('-i', '--input', required=True, help='Input FASTQ file')
    parser.add_argument('-o', '--output', required=True, help='Output FASTQ file')
    parser.add_argument('-g', '--gc', default='0,100', help='GC bounds (e.g., "30,70" or "50")')
    parser.add_argument('-l', '--length', default='0,10000', help='Length bounds (e.g., "50,100" or "200")')
    parser.add_argument('-q', '--quality', type=float, default=0, help='Minimum average quality score')
    parser.add_argument('--log', default='fastq_filter.log', help='Log file path')
    
    args = parser.parse_args()
    
    setup_logging(args.log)
    
    try:
        gc_bounds = parse_gc_bounds(args.gc)
        length_bounds = parse_length_bounds(args.length)
        
        filter_fastq(
            input_fastq=args.input,
            output_fastq=args.output,
            gc_bounds=gc_bounds,
            length_bounds=length_bounds,
            quality_threshold=args.quality
        )
    except Exception as e:
        logging.error(f"Application error: {str(e)}")

if __name__ == '__main__':
    main()

# python filter.py -i example_fastq.fastq -o output.fastq -g 30,70 -l 50,150 -q 20