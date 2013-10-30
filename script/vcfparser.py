#!/usr/bin/env python
import os, sys, csv

chrmap = {"./.":"", "0/0":"0", "0/1":"1", "1/1":"2"}

split_point = int(sys.argv[2])

#with open(sys.argv[1]) as file:
	#ks = file.next().strip().split("\t")
	#ks.remove("INFO")

with open(sys.argv[1]) as file:
	for line in file:
		if line.startswith("##"):
			continue
		elif line.startswith("#"):
			ks = line.strip().split("\t")
			print "\t".join(k for k in ks if k!="INFO")
		else:
			row = dict(zip(ks, line.strip().split("\t")))
			result = []
			for i in range(len(ks)):
				k = ks[i]
				if k=="INFO":
					continue
				elif i >= split_point:
					chrs = row[k][:3]
					result.append(chrmap[chrs])
				else:
					result.append(row[k])
			print "\t".join(result)


