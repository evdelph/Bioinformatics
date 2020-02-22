# Bioinformatics

## Project 1: Random Sequence Like Generator
### Getting Started
This code is composed of two files, a bash file "sequence_like.sh" and a python file called random_generator.py.
This script takes in an accession number and: 1) Downloads a fasta file, 2) Computes the composition of the sequence, 3) Generates a random sequence with the same length and composition, 4) Outputs the new sequence in a fasta file format.
```
./sequences_like 'YOUR ACCESSION NUMBER' 'N'
```
### Libraries and Dependencies
Below outlines the libraries used:
### Libraries
```
# import sys for passing fasta file #
import sys

# import random for randomizing sequence #
import random as r

# import bipython for fasta file parsing #
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
```
