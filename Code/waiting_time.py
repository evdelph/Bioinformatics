#! /home/i519/anaconda3/bin/python

import sys
import random as r

file = sys.argv[1]
motif = sys.argv[2]
comp_ratio = dict()

# open file #
with open(file,'r') as f:
	for line in f.readlines():
        	line = line.split()
        	line[1] = float(line[1])
        	line[0] = line[0].upper()
        	comp_ratio[line[0]] = line[1]
f.close()

# check if probabilties add up to one #
assert sum(comp_ratio.values()) == 1

# compute waiting time #
def waiting_time(composition, motif):
    	sequence = ''
    	while len(sequence) < len(motif) or sequence[-len(motif):] != motif:
        	ran = r.random()
        	threshold = 0
        	for nucleo in composition:
            		threshold += composition[nucleo]
            		if ran < threshold:
                		sequence += nucleo
                		break
    	return len(sequence)

# compute mean waiting time #
def mean_waiting_time(composition,motif,n):
    	w = [waiting_time(composition,motif) for i in range(n)]
    	return sum(w)/len(w)

# compute mean waiting time at convergence #
def converge(composition,motif,n):
	n_two = n*2
	v1 = mean_waiting_time(composition,motif,n)
	v2 = mean_waiting_time(composition,motif,n_two)
	f = lambda x,y : abs((x-y)/y) > .005
	while f(v2,v1):
		v1 = mean_waiting_time(composition,motif,n)
		v2 = mean_waiting_time(composition,motif,n_two)
	return v1

# function call #
print(converge(comp_ratio,motif,500))
