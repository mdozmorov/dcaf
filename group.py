#!/usr/bin/env python
from itertools import groupby
import sys

rows = (line.strip().split("\t") for line in sys.stdin)
for k, grp in groupby(rows, lambda r: [r[0],r[1],r[2]]):
    duplicates = list(grp)
    result = []
    for col in zip(*duplicates):
        if len(set(col)) == 1:
            result.append(col[0])
        else:
            result.append("|".join(col))
    print "\t".join(result)

# Comment out the melow code, used within python	
# A special version of groupby made by understanding the functions
rows = [line.strip().split("\t") for line in open('gwas808.catInfo.txt')]
# Test groupby
for diseases, snps in groupby(rows, lambda r: r[1]):
	listOfSnps = "\t".join(["%s" % snp[0] for snp in snps])
	print(diseases + "s: " + listOfSnps)
# Now, use groupby to output separate sets of SNPs
for diseases, snps in groupby(rows, lambda r: r[1]):
    listOfSnps = [] # For each disease, initialize empty list
    for snp in snps: listOfSnps.append(snp[0]) # Append to it corresponding SNPs selected by groupby
    with open(diseases + ".txt", "w") as out: # Open disease-specific file
        for snp in listOfSnps: # Write out each SNP
            out.write(str(snp) + "\n") # and include new line
