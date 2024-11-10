import os


def read_fastq(input_fastq: str) -> None:
    with open(input_fastq, 'r') as file:
        while True:
            name = file.readline().strip()
            if not name:
                break
            sequence = file.readline().strip()
            file.readline()
            quality = file.readline().strip()

            yield name, sequence, quality


def write_fastq(output_fastq: str, name: str, sequence: str, quality: str) -> None:
    output_dir = "Filtered"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_path = os.path.join(output_dir, output_fastq)

    with open(output_path, 'a') as file:
        file.write(f"{name}\n{sequence}\n+\n{quality}\n")
