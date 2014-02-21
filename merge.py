# Combines multiple GenomeRunner matrix files. Multiple columns and/or missing rows allowed
import pandas
import sys

table = None
for file in sys.argv[1:]:
	df = pandas.read_csv(file, index_col=0, sep="\t")
	table = df if (table is None) else pandas.merge(table, df, left_index=True, right_index=True, how="outer")

table.to_csv(sys.stdout, sep="\t", na_rep="NaN")
