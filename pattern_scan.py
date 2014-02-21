from pyfasta import Fasta
import re
import os, os.path
import sys

#f1=f['chr1'][20000000:20000100]
pattern = re.compile(sys.argv[2])
S_8x10-5_RFX5_100bp.bed
for chrom in f.iterkeys():
	seq = str(f[chrom]).upper()
	for m in re.finditer(pattern, seq):
		print "\t".join(map(str, [chrom, m.start(), m.end()]))

'''
for root, dirs, files in os.walk('/data/genomes/mm9/'):
	for f in files:
'''
