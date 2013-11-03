#!/usr/bin/env python
import sys

def col_counts(path):
	counts = []
	with open(path) as file:
		for row in file:
			cols = row.strip().split("\t")
			if not counts:
				counts = [0 for item in cols]
			for i,c in enumerate(cols):
				if c:
					counts[i] += 1
	return counts

def remove_rows(path, percent=0.1):
	cutoff = 0
	removed = 0
	with open(path) as file:
		for i,line in enumerate(file):
			fields = line.strip().split("\t")
			if not cutoff:
				cutoff = len(fields) * (1-percent)
			if len([f for f in fields if f]) > cutoff:
				print line,
			else:
				removed += 1
	sys.stderr.write("%s rows removed." % str(removed))

def remove_cols(path, percent=0.1):
	counts = col_counts(path)
	cutoff = (1-percent) * max(counts)
	# sys.stderr.write("%s cutoyff") % cutoff
	keep = [i for i,c in enumerate(counts) if c>cutoff]
	sys.stderr.write("%s cols will be removed.\n" % str(len(counts) - len(keep)))
	with open(path) as file:
		for line in file:
			fields = line.strip().split("\t")
			print "\t".join([fields[i] for i in keep])

if __name__ == "__main__":
	if sys.argv[1] == "rows":
		remove_rows(sys.argv[2])
	elif sys.argv[1] == "cols":
		remove_cols(sys.argv[2])

