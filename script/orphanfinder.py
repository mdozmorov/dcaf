#!/usr/bin/env python
from itertools import groupby
import sys

rows = (line.strip().split("\t") for line in sys.stdin)
for k, grp in groupby(rows, lambda r: r[3].split("/")[0]):
    lines = list(grp)
    if len(lines) == 1:
        print k
