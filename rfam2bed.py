import MySQLdb
import sys
import os

OUTDIR = sys.argv[1]
if not os.path.exists(OUTDIR):
    os.mkdir(OUTDIR)

conn = MySQLdb.connect(host="wren.omrf.org", db="hg19")
c = conn.cursor()
c.execute("select distinct (geneSymbol) from kgXref where rfamAcc!=''")
symbols = [r[0] for r in c]

for symbol in symbols:
    print symbol
    c.execute("""select chrom, txStart, txEnd from knownGene 
              inner join kgXref
              on kgXref.kgId=knownGene.name
              where kgXref.geneSymbol = '%s'""" % symbol)
    with open(os.path.join(OUTDIR, symbol.replace(" ", "_")), "w") as h:
        for row in c:
              h.write("\t".join(map(str, row)) + "\n")
    
