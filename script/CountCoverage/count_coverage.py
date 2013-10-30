#!/usr/bin/env python
import sys

last_chrom = None
last_end = None
count = 0
for line in sys.stdin:
    chrom, start, end = line.split("\t")[:3]
    start = int(start)
    if not end:
        end = int(start) + 1
    else:
        end = int(end)

    if chrom == last_chrom and start < last_end:
        count += (end - last_end)
    else:
        count += (end - start)

    last_chrom = chrom
    last_end = end
print "Count:", count
        
