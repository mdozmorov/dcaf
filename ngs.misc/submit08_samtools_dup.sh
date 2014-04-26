#!/bin/bash

submit1() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $2_samtools_keepdup "samtools view -b -o th_stranded_sort_rmdup/$2.pS.keepdup.bam -f 0x400 $1"
}

submit2() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $2_samtools_rmdup "samtools view -b -o th_stranded_sort_rmdup/$2.pS.rmdup.bam -F 0x400 $1" 
}

for file in `find th_stranded_sort_rmdup/ -name '*pS.rd.bam' | sort`; do
	sample=${file:23:14}
	echo $file, $sample
	submit1 $file $sample
	submit2 $file $sample
done
