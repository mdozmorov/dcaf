import fileinput, math, itertools
from collections import defaultdict

def subs(tokens, nmin, nmax):
   for i in range(len(tokens)):
      for j in range(i+nmin,min(i+nmax+1,len(tokens)+1)):
         seq = tuple(tokens[i:j])
         if not any(map(lambda t: t.isdigit(),seq)):
            if not "and" in seq:
               yield seq

tfreqs = defaultdict(int)
ngrams = defaultdict(int)
nmin = 3
nmax = 6

for line in fileinput.input():
   if line == "---":
      continue
   #tokens = line.replace("(","").replace(")","").split()
   tokens = line.split()
   for s in subs(tokens,nmin,nmax):
      ngrams[s] += 1
   for t in tokens:
      tfreqs[t] += 1

groups = itertools.groupby(sorted(ngrams,key=len,reverse=True),key=len)
for n,ks in groups:
   for k in ks:
      wb = tuple(k[:len(k)-1])
      wa = tuple(k[1:])
      if ngrams[k]==ngrams[wb]:
         del ngrams[wb]
      if ngrams[k]==ngrams[wa]:
         del ngrams[wa]

for lr,ng in sorted([(ngrams[ng], ng) for ng in ngrams],reverse=True):
   print " ".join(ng) + "\t" + str(lr)
