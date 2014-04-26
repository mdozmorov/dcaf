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

tr -d '\0' \
    | awk 'ORS=(NR%4==0)?"\n":"\0"' \
    | awk -F'\0' \
        'BEGIN {for (i=0; i<256; i++) {
            chmap[sprintf("%c",i)] = i; \
        }} \
        { split($4,ch,""); sum = 0; \
        for (i in ch) { \
            sum += chmap[ch[i]]
        }; print $1"\0"$2"\0"sum"\0"$3"\0"$4}' \
    | sort -t"\0" -r -k 2,2 -n -k3,3 \
    | sort -t"\0" -u -k 2,2 \
    | cut -d "" -f1-2,4-5 \
    | tr '\0' '\n'
