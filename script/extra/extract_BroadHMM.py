#!/usr/bin/env python2

import mysql.connector
import itertools
import operator
import sys
import os.path

OUTDIR = sys.argv[1]

class UCSCDB(object):
    def __init__(self, user, password, host, db):
        self.db = mysql.connector.connect
		
    # Function to process one table at a time
    def save_table_grouped_by_name(self, table, out_dir):
        self.c = self.db.cursor()
        self.q = """SELECT chrom, chromStart, chromEnd, name FROM %s ORDER BY name;"""
        self.c.execute(self.q, (table,))
        for grp, rows in itertools.groupby(self.c, operator.itemgetter(3)):
            path = os.path.join(OUTDIR, grp.replace("/","")+".bed")
            with open(path, "w") as h:
                for row in rows:
                    print >> h, "\t".join(map(str, row))
	
	# Function fo process alike tables
	def get_tables(self, like):
		self.c = self.db.cursor()
		self.q =  """SHOW TABLES LIKE %s;"""
		self.c.execute(self.q, (like,))
		tables = self.c.fetchall()[0]
		for table in tables:
			db.save_table_grouped_by_name(table, OUTDIR)

db = UCSCDB('genome','','genome-mysql.cse.ucsc.edu','hg19') # UCSC connection settings
# db = UCSCDB('genomerunner','genomerunner','wren2','hg19')
# db.get_tables('wgEncodeBroadHmm%')
db.save_table_grouped_by_name('wgEncodeBroadHmmGm12878HMM', OUTDIR)
    


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
    
