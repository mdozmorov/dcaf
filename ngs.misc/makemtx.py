#!/usr/bin/env python
import sys, os

folders = sys.argv[1:]

expression = {}
for folder in folders:
	with open(os.path.join(folder, "genes.fpkm_tracking")) as lines:
		header = lines.next()
		e = {}
		for line in lines:
			fields = line.strip().split("\t")
			id = fields[0]
			fpkm = fields[9]
			e[id] = float(fpkm)
		expression[folder] = e

#for folder in expression:
#	e = expression[folder]
#	for id in e:
#		transcripts.add(id)
#
#for folder,e in expression.items():

transcripts = set()
for e in expression.values():
	for id in e:
		transcripts.add(id)

print "\t".join(folders)
for tx in transcripts:
	row = [tx]
	row += [str(expression[f].get(tx,0)) for f in folders]
	print "\t".join(row)
