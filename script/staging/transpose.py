#!/usr/bin/env python

import sys, os

with open(sys.argv[1]) as file:
	outfile = sys.argv[2]
	headers = []
	vs = []
	for line in file:
		headers.append(line.split("\t")[0].strip())
		vs.append(line.split("\t")[1].strip())
	if not os.path.exists(outfile):
		with open(outfile, "w") as out:
			out.write("\t".join(headers[1:]) +"\n")
	with open(outfile, "a") as out:
		out.write("\t".join(vs) + "\n")
