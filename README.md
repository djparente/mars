# MARS - Prot
### Maintainer of Alignments using Reference Sequences for Proteins

MARS-Prot is a maintenance tool for protein multiple sequence alignments (MSAs) that allows
for the incorporating of new sequences into existing (hand-optimized) alignments.

Please see the [main documentation file](doc/MARS_readme.pdf) for information on
usage, dependencies, file formats and troubleshooting.

Directory structure:
- src - C# source code
 - mars - MARS-Prot
 - msa-comp - Tool to compare an alignment against a "gold standard"
 - DJPCommonBioinfo - Bioinformatics library that mars and msa-comp depend upon
- bin - Windows binaries
- doc - Documentation
 - doc/MARS_readme.pdf - Main documentation file
- biology - Biological data
 - LacI-GalR_Expanded.fasta - Expanded alignment of the LacI/GalR transcription regulators
 - BLOSUM62.txt - BLOSUM62 similarity matrix from NIH NCBI
