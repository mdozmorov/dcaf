#!/bin/bash
i=0
for interval in 'cat $1'; do
	i=$i+1
	for file in *.bam; do
        	samtools view -b $file $interval -o $file-$i
	done
done



