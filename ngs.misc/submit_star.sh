#!/bin/bash
submit1() {
	qsub -b y -pe threaded 4-24 -l h_vmem=8G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $1 "tophat -p \$NSLOTS -G /home/dozmorovm/work/mm9.genes.gtf --no-coverage-search --microexon-search --prefilter-multihits --b2-very-fast  --library-type=fr-unstranded -o th_unstranded/$1 /Volumes/hts_core/Shared/bowtie2_index/mm9/mm9 $2"
}

submit1 test_th-unstranded 1.gz
exit

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

