#! /home/i519/anaconda3/bin/python# 
Import libraries 
import argparse
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
from Bio.Align import PairwiseAligner
import numpy as np

# Add parse arguments
parser = argparse.ArgumentParser(description='Computes a pairwise similarity matrix from a fasta file.')
parser.add_argument('-f',help='name of the fasta file',required=True)
parser.add_argument('-s',help='name of subsitution matrix from BioPython',required=True,choices=MatrixInfo.available_matrices)
parser.add_argument('-go',help='gap opening score',type=float,required=True)
parser.add_argument('-ge',help='gap extension score',type=float,required=True)
args = parser.parse_args()

# Parse fasta file #
seqs = list(SeqIO.parse(args.f,'fasta'))

# Get substitution matrix 
substitution_matrix = getattr(MatrixInfo,args.s)

#Pairwise alignment 
aligner = PairwiseAligner()
aligner.open_gap_score, aligner.extend_gap_score  = args.go, args.ge
aligner.substitution_matrix = substitution_matrix

# Align sequences and build matrix
def similarity_matrix(seqs,n=len(seqs)):
  similarity_matrix = np.zeros([n,n])
	for i in range(len(seqs)):
    for j in range(len(seqs)):
			alignment = aligner.align(seqs[i].seq,seqs[j].seq)
			similarity_matrix[i][j] = alignment.score
	return similarity_matrix

m = similarity_matrix(seqs)

def print_matrix(m):
	for i in m:
		row = ""
		for j in i:
			row+=str(j)+'\t'
		  print(row)

print_matrix(m)
