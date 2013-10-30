import fileinput, collections
from contextlib import closing

d = collections.defaultdict(int)
with closing(fileinput.input()) as lines:
	for line in lines:
		fields = line.strip().split("\t")
		d[fields[5]]+= 1
		print line.strip()+"_"+str(d[fields[5]])
		
		
