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

###### $\color{red}{\text{NEW !}}$ `filter_fastq`
- **Filtering sequences** by:
   - GC-content
   - average reading quality
   - length

<div align="center"> <h2 align="center"> Installation </h2> </div>

Clone the repository: `git clone git@github.com:GitArthurLee/bioinflib.git`

<div align="center"> <h6 align="center"> Example of using </h6> </div>

``` python
from bioinflib import filter_fastq

EXAMPLE_FASTQ = {'name' : ('sequence', 'quality')}
print(filter_fastq(EXAMPLE_SEQS, gc_bounds = 60, ength_bounds = 2**32, quality_threshold = 0))
```
<div align="center"> <h2 align="center"> Contact </h2> </div>

Created by aal1999arth@gmail.com

Institute of Bioinformatics
