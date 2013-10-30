#!/usr/bin/env python
import csv
import sys
import os
import shutil

# a=[]
# reader = csv.reader(open(sys.argv[1],"r"), delimiter='\t')
# for row in reader:
# 	a.append(row)
# print a

with open("gwas2bed_02-26-2013_categories.txt") as h:
	for line in h:
		file, type = line.split("\t")
		dest = os.path.join(type.strip().lower(), file)
		src = os.path.join("gwas2bed_02-26-2013", file)
		print src
		if os.path.exists(src):
			shutil.copyfile(src, dest)
		

#data_files = [(x[0], x[2]) for x in os.walk('gwas2bed_02-26-2013')]
#print data_files[0][1]


