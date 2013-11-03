#!/usr/bin/env python
from itertools import groupby
import sys

rows = (line.strip().split("\t") for line in sys.stdin)
for k, grp in groupby(rows, lambda r: [r[0],r[1],r[2]]):
    duplicates = list(grp)
    result = []
    for col in zip(*duplicates):
        if len(set(col)) == 1:
            result.append(col[0])
        else:
            result.append("|".join(col))
    print "\t".join(result)
