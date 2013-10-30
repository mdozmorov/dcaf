#!/usr/bin/env python
import sys
from operator import itemgetter

with open(sys.argv[1]) as lines:
	lines = lines.readlines()
	hdr = [line.split()[1:] for line in lines[:4]]
	names = ["\t".join(r) for r in zip(*hdr)]
	for line in lines[5:]:
		fields = line.split("\t")
		feature = fields[0]
		print feature
		counts = map(int, fields[1:])
		with open(sys.argv[2]+"/"+feature+".bed", "w") as out:
			nc = sorted([(n,c) for (n,c) in zip(names,counts) if c>0], key=itemgetter(1))
			out.write("\n".join([n+"\t"+`c` for (n,c) in reversed(nc)]))
