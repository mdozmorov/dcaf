# Reads 3-column file containing "line #", "Gene(s)", "Ratio"
# Genes may be " /// " separated
# Splits gene line into a set of genes
# Outputs "Gene" - "Ratio" pair

import fileinput

h = fileinput.input()
h.__next__()

for line in h:
    fields = line.replace('"', '').strip().split('\t')
    fc = fields[2]
    genes = fields[1].split(" /// ")
    
    for gene in genes:
        print(gene, fc, sep="\t")
