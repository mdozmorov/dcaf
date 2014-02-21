import sys
for line in sys.stdin:
	if line.startswith("##") or not line.strip():
		print line,
	else:
		fields = line.strip().split("\t")
		fields[2] = fields[2].split(";")[0]
		gene_id = fields[-1].split(";")[0]
		fields[-1] = "gene_id \"%s\"; transcript_id \"%s.1\";" % (gene_id,gene_id)
		tmp = fields[7]
		fields[7] = fields[5]
		fields[5] = tmp
		print "\t".join(fields)
