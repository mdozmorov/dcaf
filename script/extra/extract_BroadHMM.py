#!/usr/bin/env python2

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
    def save_table_grouped_by_name(self, table, out_dir):
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        c = self.db.cursor()
        q = """SELECT chrom, chromStart, chromEnd, name FROM %s ORDER BY name;""" % table
        c.execute(q)
        for grp, rows in itertools.groupby(c, operator.itemgetter(3)):
            path = os.path.join(out_dir, grp.replace("/","")+".bed")
            with open(path, "w") as h:
                for row in rows:
                    print(*row, sep="\t", file=h)

    # Function fo process alike tables
    def get_tables(self, like, out_dir):
        c = self.db.cursor()
        q =  """SHOW TABLES LIKE %s;"""
        c.execute(self.q, (like,))
        tables = c.fetchall()[0]
        for table in tables:
            db.save_table_grouped_by_name(table, out_dir)

def main():
    OUTDIR = sys.argv[1]
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

if __name__ == "__main__":
    main()


