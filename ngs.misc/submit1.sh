#!/bin/bash

submit1() {
	qsub -b y -pe threaded 4-24 -l h_vmem=8G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $1 "tophat -p \$NSLOTS --mate-inner-dist 200 --mate-std-dev 50 --no-coverage-search --microexon-search --prefilter-multihits --b2-very-sensitive -z gzip --library-type=fr-unstranded --transcriptome-index=/home/dozmorovm/work/transcriptome_data/known -o $1 /Volumes/hts_core/Shared/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome $2 $3"
}

k=0
for file in `ls *.gz | sort`; do
	if [ $k = 0 ]
	then
		s1=$file
		k+=1
	else
		s2=$file
		k=0
		sample=${s1%????????????????????????????}
		submit1 $sample $s1 $s2
		sleep 5
	fi
done
