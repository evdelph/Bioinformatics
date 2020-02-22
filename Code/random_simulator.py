#! /home/i519/anaconda3/bin/python

# import sys for passing fasta file #
import sys

# import random for randomizing sequence #
import random as r

# import bipython for fasta file parsing #
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

# assign input variables from bash script #
fasta_file = sys.argv[1]
n = int(sys.argv[2])

# parse fasta file #
seq_file = list(SeqIO.parse(fasta_file,"fasta"))

# calculate composition and ratio #
alphabet = 'ATCG'
seq = [y for x in seq_file[0].seq for y in x if y in alphabet]
seq_dict = dict(map(lambda nucleo : (nucleo,seq.count(nucleo)),set(seq)))
seq_ratio =  dict(map(lambda nucleo : (nucleo,round(seq.count(nucleo)/len(seq),2)),set(seq)))

assert abs(sum(seq_dict.values()) - len(seq_file[0].seq)) < .05
assert abs(sum(seq_ratio.values()) - 1) < .05

# generate n random sequences #
for i in range(n):
	seq = list(seq_file[0].seq)
	# generate random sequence #
	r.shuffle(seq)
	seq_file[0].seq  = "".join(seq)
	# create new record with new random sequence #
	record = SeqRecord(Seq(seq_file[0].seq,IUPAC.unambiguous_dna),seq_file[0].id,seq_file[0].name,seq_file[0].description)
	print(record.format("fasta"))
