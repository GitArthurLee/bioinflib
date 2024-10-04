<div align="center"> <h1 align="center"> $\color{green}{\text{BIOINFLIB}}$ </h1> </div>

<div align="center"> <h6 align="center"> Arthur's first bioinformatician library </h6> </div>

<div align="center"> <h3 align="center"> Tools for analyzing and manipulating DNA and RNA sequences, including transcription, complementation, and sequence filtering. </h3> </div> 

<div align="center"> <h2 align="center"> Features </h2> </div>

- DNA transcription
- RNA transcription
- Creating complementary sequences
- Reversing sequences
- Calculation of the molecular weight of single-stranded DNA and RNA
- Counting the GC content
- Calculation of the melting point of primers (Tm)

###### NEW
- **Filtering sequences** by:
   - GC content
   - average reading quality
   - length

<div align="center"> <h2 align="center"> Installation </h2> </div>

Clone the repository: `git clone git@github.com:GitArthurLee/bioinflib.git`

###### Example
``` python
from bioinflib import filter_fastq

EXAMPLE_FASTQ = {'name' : ('sequence', 'quality')}
print(filter_fastq(EXAMPLE_SEQS, gc_bounds = 60, ength_bounds = 2**32, quality_threshold = 0))
```
<div align="center"> <h2 align="center"> Contact </h2> </div>

Created by aal1999arth@gmail.com

