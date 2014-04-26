#!/bin/bash

# Usage:
# ./btang.sh < my.fq > out.fq
# zcat my.fq.gz | ./btang.sh > out.fq

# Removes duplicate sequences from FASTQ. How it works:
# * Convert each FQ record to a one-per-line format, 
#   seperated by null characters
# * Create a temporary field, 'sum', the sum of the quality scores
# * Sort, first by sequence then by quality sum
# * Take the first unique sequence (which will have the highest quality)
# * Convert back to regular FASTQ

zcat $1 \
	| tr -d '\0' \
	| awk 'ORS=(NR%4==0)?"\n":"\0"' \
    	| sort -k1,1 > $$tmp1
zcat $2 \
	| tr -d '\0' \
	| awk 'ORS=(NR%4==0)?"\n":"\0"' \
    	| sort -k1,1 > $$tmp2
join -t '\0' -1 1 -2 1 $$tmp1 $$tmp2 > $$joined
cut -d "" -f1-4 | tr '\0' '\n' > `basename $1 .fastq.gz`.sorted.fastq
cut -d "" -f5-8 | tr '\0' '\n' > `basename $1 .fastq.gz`.sorted.fastq
