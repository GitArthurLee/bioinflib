<div align="center"> <h1 align="center"> $\color{green}{\text{BIOINFLIB}}$ </h1> </div>

<div align="center"> <h6 align="center"> Arthur's first bioinformatician library </h6> </div>

<div align="center"> <h3 align="center"> Tools for analyzing and manipulating DNA and RNA sequences, including transcription, complementation, and sequence filtering. </h3> </div> 

<div align="center"> <h2 align="center"> Features </h2> </div>

###### `run_dna_rna_tools`
- DNA and RNA transcription `trancribe`
- Creating complementary sequences `complement`
- Reversing sequences `reverse`
- Calculation of the molecular weight of single-stranded DNA and RNA `MW`
- Counting the GC-content `GC`
- Calculation of the melting point of primers for PCR `tm`

###### `filter_fastq` 🔥 $\color{pink}{\text{Now with .FASTQ reading support}}$ 🔥
- **Filtering sequences** by:
   - GC-content
   - average reading quality
   - length

###### $\color{red}{\text{NEW!}}$ `convert_multiline_fasta_to_oneline`
- Сonverts a multi-line sequence into a single line

###### $\color{red}{\text{NEW!}}$ `parse_blast_output`
- Parses the BLAST output, extracting a description of the first match for each request.

<div align="center"> <h2 align="center"> Installation </h2> </div>

Clone the repository: `git clone git@github.com:GitArthurLee/bioinflib.git`

<div align="center"> <h6 align="center"> Example of using </h6> </div>

``` python
from bioinflib import filter_fastq

print(filter_fastq('example.fastq', gc_bounds = 60, ength_bounds = 2**32, quality_threshold = 0))
```

``` python
import bio_files_processor 

convert_multiline_fasta_to_oneline('example_multiline.fasta')
parse_blast_output('example_blast_results.txt')
```
<div align="center"> <h2 align="center"> Contact </h2> </div>

Created by aal1999arth@gmail.com 💚

Bioinformatics Institute, 2024