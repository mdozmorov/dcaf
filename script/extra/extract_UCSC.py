#!/usr/bin/env python2
## ==================================================================
## Get data from UCSC MySQL server, extract hg19 genomic coordinates in separate
## files defined by 'field' column
##
## Usage: python extract_UCSC.py [table name] [field name] [output dir]
##
## Example: python extract_UCSC.py wgEncodeBroadHmmGm12878HMM histoneMarks
## Example: python extract_UCSC.py wgEncodeRegTfbsClusteredV3 name tfbsEncode
## Example: python extract_UCSC.py gap type gapLocations
## Example: python extract_UCSC.py coriellDelDup cellType coriellCNVs
## Example: python extract_UCSC.py gwasCatalog trait gwasCatalog
## Example: python extract_UCSC.py knownAlt name altSplicing
## Example: python extract_UCSC.py wgRna type ncRnas
## Example: python extract_UCSC.py tfbsConsSites name tfbsConserved
## Example: python extract_UCSC.py dgvMerged varType genomicVariants
## Example: python extract_UCSC.py nestedRepeats repClass repeats
## Example: python extract_UCSC.py snp138Flagged func snpClinical
## =================================================================
import mysql.connector
import itertools
import operator
import sys
import os.path

class UCSCDB(object):
    def __init__(self, user, password, host, db):
        self.db = mysql.connector.connect(user=user, 
                password=password, host=host, database=db)

    # Function to process one table at a time
    def save_table_grouped_by_name(self, table, tname, out_dir):
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        c = self.db.cursor()
        q = """SELECT chrom, chromStart, chromEnd, name, %s FROM %s ORDER BY %s;""" % (tname, table, tname)
        c.execute(q)
        for grp, rows in itertools.groupby(c, operator.itemgetter(4)):
            path = os.path.join(out_dir, grp.replace("/","")+".bed")
            with open(path, "w") as h:
                for row in rows:
                    print(*row[:4], sep="\t", file=h)

    # Function fo process alike tables
    def get_tables(self, like, out_dir):
        c = self.db.cursor()
        q =  """SHOW TABLES LIKE %s;"""
        c.execute(self.q, (like,))
        tables = c.fetchall()[0]
        for table in tables:
            db.save_table_grouped_by_name(table, out_dir)

def main():
    TABLE = sys.argv[1]
    TNAME = sys.argv[2]
    OUTDIR = sys.argv[3]
    db = UCSCDB('genome','','genome-mysql.cse.ucsc.edu','hg19') # UCSC connection settings
    # db = UCSCDB('genomerunner','genomerunner','wren2','hg19')
    # db.get_tables('wgEncodeBroadHmm%')
    db.save_table_grouped_by_name(TABLE, TNAME, OUTDIR)
        


    # q = """SHOW TABLES LIKE wgEncodeBroadHmm%"""  

    ## Works for a single table
    # db = mysql.connector.connect(user="genomerunner", 
        # password="genomerunner", database="hg19", host="wren2")
    # c = db.cursor()

    # q = """SELECT chrom, chromStart, chromEnd, name FROM wgEncodeBroadHmmGm12878HMM ORDER BY name;"""
    # c.execute(q)

    # for grp, rows in itertools.groupby(c, operator.itemgetter(3)):
        # path = os.path.join(OUTDIR, grp.replace("/","")+".bed")
        # with open(path, "w") as h:
            # for row in rows:
                # print >> h, "\t".join(map(str, row))

if __name__ == "__main__":
    main()


