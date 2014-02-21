#!/usr/bin/env python
import sys, re
for path in sys.stdin:
	txt = open(path.strip()).read()	
	seqs = re.findall("Overrepresented(.+?)>>END", txt, re.DOTALL)
	if not seqs:
		continue
	path_out = "LB_" + re.findall("LB_(.+?)_fastqc", path)[0]
	for line in seqs[0].split("\n")[2:]:
		print line.strip() + "\t" + path_out

