#!/usr/bin/env python

import sys

PREFIX = sys.argv[1]

with open(PREFIX+".map") as lines:
	columns = zip(*[line.strip().split() for line in lines])
	print "\t".join(columns[1])	

with open(PREFIX+".ped") as lines:
	for line in lines:
		fields = line.strip().split()
		out = [fields[0]]
		for i in range(6, len(fields), 2):
			al1, al2 = fields[i:i+2]
			if (al1 == "0") or (al2 == "0"):
				out.append("0")
			elif al1 == al2:
				out.append("2")
			else:
				out.append("1")
		print "\t".join(out)
