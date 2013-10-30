#!/usr/bin/env python
# Parse .BED records by length, output in different files

from itertools import groupby
import fileinput
# Define length of genomic interval
breaks = [24, 100, 200, 500, 1000, 10000]

def read_line_categories(lines):
	for line in lines:
		fields = line.split("\t")
		length = int(fields[2]) - int(fields[1])
		category = len(breaks)
		for i, n in enumerate(breaks):
			if length < n:
				category = i
				break
		yield (category, line)

line_categories = read_line_categories(fileinput.input())
for cat, lines in groupby(line_categories, lambda line: line[0]):
	with open(str(cat)+".txt", "a") as out:
		for _,line in lines:
			out.write(line)
