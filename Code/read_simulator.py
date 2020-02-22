#! /home/i519/anaconda3/bin/python

# Import libraries #
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Seq import MutableSeq 
import random as rn
import math

parser = argparse.ArgumentParser(description='Read simulator.')
parser.add_argument('-g',help='reference genome',required=True)
parser.add_argument('-m',help='mutation',type=float,required=True)
parser.add_argument('-n',help='number of reads',type=float,required=True)
parser.add_argument('-l',help='length of reads',type=float,required=True)
parser.add_argument('-e',help='error rate',type=float,required=True)
parser.add_argument('-o',help='output (fasta file)',required=True)
args = parser.parse_args()


# Parse fasta file #
input = list(SeqIO.parse(args.g,'fasta'))
seqs = MutableSeq(str(input[0].seq[:1000000]))

# Function to introduce mutations # 
def mutate_genome(seq,m):
	mutation_rate = int(m * len(seq))
	for i in range(mutation_rate):
		c = rn.randint(0,len(seq)-1)
		seq[c] = mutate(seq[c])
	return seq

# Pick a character to mutate #
def mutate(char):
         mutation_set = list({'A','G','C','T'} - {char})
         return mutation_set[rn.randint(0,len(mutation_set)-1)]

# Introduce errors and write reads to new file #
def create_reads(seq,n,l,e,o):
    # Round up error rates #
	error_rate = int(math.ceil(e * l))
	with open(o, "w") as output_handle:
		for i in range(n):
			# Get sample read #
			index = rn.randint(0,len(seq)-1)
			sample_read  = ""
			if index + l > len(seq):
				start = index-((index + l) - len(seqs))
				sample_read = seq[start:start+l]
			else:
				sample_read = seq[index:index+l]
			# Introduce technical errors #
			for y in range(error_rate):
				c = rn.randint(0,l-1)
				hold = sample_read[c]
				sample_read[c] = mutate(sample_read[c])
			# Read mutated read to ouput file #
			read = SeqRecord(Seq(str(sample_read)), name = "Read: "+str(i+1),id = str(i+1), description = "Sample read "+str(i+1) + " out of "+str(n))
			SeqIO.write(read, output_handle, "fasta")


mutate_seq  = mutate_genome(seqs,args.m)
create_reads(mutate_seq,int(args.n),int(args.l),args.e,args.o)
