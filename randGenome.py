#!/usr/bin/python
## ==================================================================
## Generates random sequence of letters from the list. Length is the first argument
##
## Usage: python randGenome.py 10 > randomGenome10.txt
## ==================================================================
import random
import sys

SEQ=['A', 'T', 'C', 'G']
for i in range(int(sys.argv[1])):
	sys.stdout.write(SEQ[random.choice([0,1,2,3])])