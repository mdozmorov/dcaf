from collections import defaultdict

MM = {}
with open("GSE15907new-mtx.txt") as h:
    header = h.next() 
    for line in h:
        fields = line.strip().split("\t")
        MM[fields[0]] = map(float, fields[1:])

human2mouse = defaultdict(set)
with open("Cross-species - HMD_Human5.rpt") as h:
    ok = False
    for line in h:
        if not ok and line.startswith("Human Chr"):
            ok = True
        elif ok:
            fields = line.strip().split("\t")
            if not len(fields) == 10:
                continue
            human = fields[2]
            mouse = fields[7]
            human2mouse[human.upper()].add(mouse.upper())

symbol2probe = defaultdict(set)

with open("GPL6246-21513.txt") as h:
    ok = False
    for line in h:
        fields = line.strip().split("\t")
        if line.startswith("ID"):
            ok = True
        elif ok and len(fields) == 12:
            if fields[9] == "---":
                continue
            affy = fields[0]
            symbol = fields[9].split("//")[1].strip()
            symbol2probe[symbol.upper()].add(affy)

HS = {}
with open("ShannonData1.txt") as h:
    hdr = h.next()
    for line in h:
        fields = line.strip().split("\t")
        id = fields[0].replace('"',"")
        if id.upper() in human2mouse:
            HS[id.upper()] = map(float, fields[1:])

human2affy = defaultdict(set)
for hg, mms in human2mouse.items():
    for mm in mms:
        if mm in symbol2probe:
            print hg, mm, symbol2probe[mm]
        for affy in symbol2probe[mm]:
            human2affy[hg].add(affy)

print len(human2affy)

h = open("matrix.tab", "w")
for human_gene in HS:
    hs_expr = HS[human_gene]
    if not human_gene in human2affy:
        continue
    mouse_expr = sorted([MM[affy] for affy in human2affy[human_gene]], key=sum, reverse=True)[0]
    expr = [human_gene] + hs_expr + mouse_expr
    print >> h, "\t".join(map(str, expr))

h.close()    
