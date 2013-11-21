#!/usr/bin/env python2

import sys

with open(sys.argv[1]) as file:
	for line in file:
		fields = line.strip().split("\t")
		chrom, start, end = fields[:3]
		start = int(start)
		end = int(end)
		if ((end - start) == 1):
			print "\t".join([chrom, str(start), str(end)])
		else:
			print "\t".join([chrom, str(start), str(start+1)])
			