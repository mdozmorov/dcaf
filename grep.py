import sys
import re

def grep(lines, patterns):
    return filter(lambda line: any(p.match(line) for p in patterns), lines)

    #return (line for line in lines
    #        if any(p.match(line) for p in patterns))

    #for line in lines:
    #    if any(p.match(line) for p in patterns)):
    #        yield line

patterns = list(map(re.compile, sys.argv[1:]))

for line in grep(sys.stdin, patterns):
    print(line)
