#!/bin/bash
submit1() {
	qsub -b y -pe threaded 4-24 -l h_vmem=8G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $1 "tophat -p \$NSLOTS -G /home/dozmorovm/work/mm9.genes.gtf --mate-inner-dist 200 --mate-std-dev 50 --coverage-search --microexon-search --prefilter-multihits --b2-very-sensitive -z gzip --library-type=fr-secondstrand 	-o th_stranded/$1 /Volumes/hts_core/Shared/bowtie2_index/mm9/mm9 $2 $3"
}

#submit1 test_th-stranded 1.gz
#exit

k=0
for file in `find . -name *.gz -type f | sort`; do
	if [ $k = 0 ]
	then
		s1=$file
		k+=1
	else
		s2=$file
		k=0
		sample=`basename $s1`
		sample=${sample:0:17}
		echo $sample $s1 $s2
		submit1 $sample $s1 $s2
		sleep 1
	fi
done

