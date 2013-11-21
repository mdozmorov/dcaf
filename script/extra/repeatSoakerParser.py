#!/usr/bin/python
import gzip
from random import choice
from collections import namedtuple
import sys
## ==================================================================
## Converts gzipped FASTA file with sequences into contiguous sequence,
## interspersed with random sequences. Output into stdout
##
## Usage: python repeatSoaker.py [filename.fa.gz] > [filenameOut.fa]
## ==================================================================

FASTA = namedtuple("FASTA", "key seq") # Generator to keep 
# SEQ = ['A', 'T', 'C', 'G']
SEQ = ['a', 't', 'c', 'g']

def read_fasta(handle):
	key = None
	for line in handle:
		if line.startswith(">"):
			if key:
				yield FASTA(key, seq)
			key = line.strip()[1:]
			seq = ""
		else:
			seq += line.strip()
			
def random_seq(n = 20): # Random sequence generator
	random_seq=""
	for i in range(0, n):
		random_seq += choice(SEQ)
	return random_seq

print(">",sys.argv[1]) # Header is the name of original file
with gzip.open(sys.argv[1], 'rt') as h:
	for record in read_fasta(h):
		print(random_seq(), end="") # Start with random sequence
		print(record.seq, end="") # Stitch sequences together
print(random_seq()) # End with random sequence