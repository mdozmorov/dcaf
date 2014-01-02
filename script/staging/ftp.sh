#!/bin/sh

#lftp ftp-trace.ncbi.nih.gov/1000genomes/data/ftp/ -e find > paths

for accession in `grep "\.mapped\." paths | grep exome_alignment | grep bam$ | grep -F ASW.txt`; do
	echo $accession
done
#	URL=`ftp ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/$accession/exome_alignment/`
#done
