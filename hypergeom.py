'''
Analysis of several FOI files against several GFs using Fisher's exact test. Best used for SNP set analysis, using whole SNP database as a spot background.
sys.argv[1] - text file with FOI file names. Include full, or relative path, if needed.
sys.argv[2] - text file with GF file names. Include full, or relative path, if needed.
sys.argv[3] - spot background file
'''

import sys
import gzip
import pybedtools
import multiprocessing as mp
import math
from scipy.stats import hypergeom
import scipy
import pdb
import numpy

# Collect lines from a file into a dictionary
def read_lines(path):
	elems = []
	with open(path) as h:
		for line in h:
			elems.append(line.strip())
	return elems

# Run hypergeometric test (pyGenome, pyFOI, pyGF)
def run_hypergeometric(g, f, b):
	foi = pybedtools.BedTool(f)
	gf = pybedtools.BedTool(g)
	bg = pybedtools.BedTool(b)
	foi_obs = len(foi.intersect(gf, u = True)) # number of FOIs overlapping with a GF
	bg_obs = len(bg.intersect(gf, u = True)) # number of spot bkg overlapping with a GF
	rnd_obs = (len(foi)*bg_obs/len(bg)) # Mean of hypergeometric distiribution
	# pdb.set_trace()
	if foi_obs == rnd_obs: # No difference
		return 0
	elif foi_obs < rnd_obs: # Underrepresentation
		# pdb.set_trace()
		# pval = hypergeom.cdf(foi_obs,len(bg),bg_obs,len(foi)) 
		# if pval <= 0:
		# 	pval = float(numpy.finfo(numpy.float64).tiny) # If pval is 0, set to min, to avoid log10 error 
		# elif pval >= 1:
		#	pval = 1 # Sometimes there may be an overflow in another direction
		pval = scipy.stats.fisher_exact([[foi_obs,bg_obs],[len(foi)-foi_obs,len(bg)-bg_obs]],alternative='less')[1]
		return math.log10(pval) # Calculating cdf hypergeometric distribution, not adding "-" to signify underrepresentation
	elif foi_obs > rnd_obs:	# Overrepresentation
		# pdb.set_trace()
		# pval = hypergeom.sf(foi_obs,len(bg),bg_obs,len(foi))
		# if pval <= 0:
		# 	pval = float(numpy.finfo(numpy.float64).tiny) # If pval is 0, set to min, to avoid log10 error 
		# elif pval >= 1:
		#	pval = 1 # Sometimes there may be an overflow in another direction
		pval = scipy.stats.fisher_exact([[foi_obs,bg_obs],[len(foi)-foi_obs,len(bg)-bg_obs]],alternative='greater')[1]
		return -math.log10(pval) # Calculating sf hypergeometric distribution
	else: # No difference
		return 0
		# print "e.obs : %d, rand.obs :			%d, len(back) : %d, back.obs : %d, len(e.A) : %d, hypergeomp_value : %r" % (eobs, randobs, len(back), backobs, len(eA), hypergeomp_value)

# Run hypergeometric for multiple FOIs from file		
def run_multiple_hypergeometric((gf, fois, bg)):
	result = [run_hypergeometric(gf,foi,bg) for foi in fois]
	result.insert(0, gf) # Add GF name to 0 position, the rest are p-values
	return result

if __name__ == "__main__":
	fois = read_lines(sys.argv[1]) # Text file with FOI file names 
	gfs = read_lines(sys.argv[2]) # Text file with GF file names, gzipped
	bg_path = sys.argv[3] # Path to spot background file

	pool = mp.Pool(mp.cpu_count())
	jobs = [(gf, fois, bg_path) for gf in gfs]

	print "\t".join(fois)
	for result in pool.imap(run_multiple_hypergeometric, jobs): # pool.imap. Use 'map' to use single thread
		print "\t".join(map(str, result))

