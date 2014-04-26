#!/bin/bash

submit1() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N s03_sicer_$2 "SICER.sh 02_bt2_aligned/ $2 InPut_D.bed 03_sicer/H3K4me1.fragment320/ hg19 1 200 320 0.74 600 0.01" 
}

submit1() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N s03_sicer_$2 "SICER.sh 02_bt2_aligned/ $1 InPut_D.bed 03_sicer/ hg19 1 200 150 0.74 600 0.01" 
}


for file in `ls 02_bt2_aligned/H3K4me1.bed | sort`; do
#for file in `ls 02_bt2_aligned/PU*.bam | sort`; do
        sample=`basename $file`
		sample=${sample:0:7}
        #submit1 $file $sample
		submit1 $file $sample
		#echo $file $sample
		#break
		#sleep 10
done

