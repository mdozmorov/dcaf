import MySQLdb
import sys
import os

OUTDIR = sys.argv[1]
if not os.path.exists(OUTDIR):
    os.mkdir(OUTDIR)

conn = MySQLdb.connect(host="localhost", db="hg19")
c = conn.cursor()
c.execute("select distinct (trait) from gwasCatalog")
traits = [r[0] for r in c]

for trait in traits:
	# trait.replace("'","\\'")
	print trait
	c.execute("select chrom, chromStart, chromEnd, name  from gwasCatalog where trait = %s group by chrom, chromStart, chromEnd, name", [trait])
	with open(os.path.join(OUTDIR, trait.replace("'","").replace(" ","_").replace("/","-")), "w") as h:
		for row in c:
			h.write("\t".join(map(str, row)) + "\n")
    
