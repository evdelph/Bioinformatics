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
## Project 2: Waiting Time Simulator
The file waiting_time.py takes a composition and a motif as input and prints the expected waiting time for that motif, given the composition. Waiting time is calculated until it's converged, less than or equal to 1%.

### Getting Started
First, the code will read in a text file of the alphabet of your sequence along with their ratios. An example is this:
```
A 0.3
B 0.3
C 0.4
```
Next, you call this file as such:
```
./waiting_time.py <composition_filename> motif
```
The code will keep running simulations of the experiment until estimate converges. This code also includes error checking to ensure all ratios add up to 1.

## Project 3: Pairwise Similarity Matrix
This script reads n protein sequences from a fasta file, does pairwise alignment between each possible pair, and generates a similarity matrix. Various alignment modules from Biopython (aka not homegrown/from scratch!) were used to optimize the alignments' accuracies.

### Getting started
To use the script, use the following snippet of code:
```
similarity_matrix.py -f <fasta_filename> -s <substitution matrix name> -go <gap opening score> -ge <gap extension score>
```
The code outputs the similarity matrix.

### Libraries Used
```
import argparse
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
from Bio.Align import PairwiseAligner
import numpy as np
```

## Project 4: Hidden Markov Simulator
This is a naive implementation of a hidden markov model, given two hidden and two output states. The code is found in HMM_generator.py. Sample output is shown below:
```
Iteration: 1

 Initial Probabilities: {0.72: 'E', 0.28: 'M'}
 Tau State E: {0.5: 'M'}
 Tau State M: {0.16: 'E', 0.84: 'M'}
 Output Probabilities: {0.41: 'T', 0.59: 'H'}
 Random Generated Sequence Length: 16

State Results: EMMMMMMMMMMMMMMM
{'%E': 0.06, '%M': 0.94}

Output Results: TTHTHHHHHHTTTHTH
{'%H': 0.56, '%T': 0.44}

 Iteration: 2

 Initial Probabilities: {0.27: 'E', 0.73: 'M'}
 Tau State E: {0.95: 'E', 0.05: 'M'}
 Tau State M: {0.72: 'E', 0.28: 'M'}
 Output Probabilities: {0.77: 'T', 0.23: 'H'}
 Random Generated Sequence Length: 54

State Results: MEEEEEEEEEEEEEEEEEEMEMEEEEEEEEEEEEEEEEMEEEEEEEEEEEEEEM
{'%E': 0.91, '%M': 0.09}

Output Results: TTHTTHTTTHHTTTTTTTTTHTHTTTTTTTTTTTTTHTHTHTHTTTTTTTTTTT
{'%H': 0.19, '%T': 0.81}
```
## Project 5: Mapping Reads Back to Genomes

### Part 1: Build a Read Simulator
The file read_simulator.py simulates short reads from a genome. It takes a genome (in fasta format) as an input, and then does the following:

* Introduces mutations (at the given rate )
* Sample reads of given length and number
* Introduces technical errors (at the given rate)

Finally, the reads are written to a fasta file. You call the script like this:
```
simulate_reads -g <genome_file> -m <mutation_rate> -e <error_rate> -n <number_of_reads> -l <length_of_reads> -o <output_file>
```
Biopython libraries were used for the simulation.

### Part 2: Map Simulated Reads Back to Genome
By using SAMtools/Tview, the next part was to look at the alignments. The report is found in HW4.pdf. It analyzes the mutations and technical alignments with the reads that were broken up in the previous step (part 1) were mapped back to that same genome. The report includes tview and coverage plot visualizations.

### Part 3: Assemble Genomes from Short Reads
This required using a super computer to submit a job. The results are located in the N50.fa and scafolds-1.fasta files.
