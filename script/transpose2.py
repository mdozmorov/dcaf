import sys
m = []
with open(sys.argv[1]) as file:
	for line in file:
		m.append(line.strip().split("\t"))
	for row in zip(*m):
		print "\t".join(row)
