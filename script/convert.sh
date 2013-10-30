#!/bin/bash

for file in GWAS_disease/*; do
	sample=`basename "$file"`
	sample=${sample%????}
	echo $sample.bed
	fgrep -w -f "$file" snp135_gwas.bed > "$sample".bed
done
