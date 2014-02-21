#!/usr/bin/env python2

import sys

with open(sys.argv[1]) as h:
    patterns = sorted([p.strip() for p in h])

for line in sys.stdin:
    field = line.split()[3].split("_")[0]
    c = cmp(field, patterns[0])
    #print field,c 
    if c == 1:
        patterns.pop(0)
        if not patterns:
            break
    elif c == 0:
        print line,
