#!/usr/bin/env python2
# Made after understanding exome2.py
# Parses rows of matrix and outputs column values meeting a criterium

import csv, sys,  argparse, os

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
        
        for row in h:
            cells = row.split()
            id = cells[0]
            if args.verbose:
                print "* Processing patient %s ..." % id
            bed_path = os.path.join(args.outdir, id+".txt")
            with open(bed_path, "w") as out:
                for snp, value in zip(header[1:], cells[1:]):                    
					if value not in ('0', '1', '2', '3', '4'):
						continue
					if value != '2':
						out.write(snp+"\n")