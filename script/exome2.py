#!/usr/bin/env python2

import csv, sys, mysql.connector, argparse, os

db = mysql.connector.connect(user="genomerunner", 
    password="genomerunner", database="hg19", host="wren2")
c = db.cursor()

def get_coordinates(snp, verbose=False):
    if snp.startswith("rs"): # rd IDs can be converted
        id = snp.split("_")[0] 
        q = """SELECT chrom, chromStart, chromEnd from snp135 WHERE name=%s"""
        c.execute(q, (id,))
        try:
            chr, start, end = c.fetchall()[0]
        except:
            if verbose:
                print >> sys.stderr, "%s" % snp # No match found SNPs
            return None
        return chr, start, end
    else:
        base = snp.split("_")[0]
        chr, start = base.split(":")
        return ''.join(['chr',str(chr)]), int(start), int(start)+1

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", "-o", default="bed")
    parser.add_argument("snp_file", nargs=1)
    parser.add_argument("--verbose", "-v", default=False,
        action="store_true")
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    with open(args.snp_file[0]) as h:
        header = list(h.next().strip().split())
        
        coordinates = dict((snp,get_coordinates(snp, verbose=args.verbose))\
             for snp in header[6:])
        for row in h:
            cells = row.split()
            id = cells[1]
            if args.verbose:
                print "* Processing patient %s ..." % id
            bed_path = os.path.join(args.outdir, id+".bed")
            with open(bed_path, "w") as out:
                for snp, value in zip(header[6:], cells[6:]):
                    if value not in ('1', '2'):
                        continue
                    if not coordinates[snp]:
                        continue
                    chr, start, end = coordinates[snp]
                    fields = [chr,start,end,snp,value,"."]
                    out.write("\t".join(map(str, fields))+"\n")

db.close()
