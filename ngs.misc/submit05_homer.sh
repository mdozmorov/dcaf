#!/bin/bash

submit1() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N s03_homer_$2 "makeTagDirectory 03_homer/tag_$2 $1 " 
}

submit2() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N s03_homer2_$2 "makeTagDirectory 03_homer/tag_$2 -genome hg19 -checkGC $1 " 
}

submit3() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N s03_homer3_$2 "makeUCSCfile 03_homer/tag_$2 -o auto" 
}

submit3.1() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N s03_homer3.1_$2 "makeUCSCfile 03_homer/tag_$2 -i 03_homer/tag_InPut_D -o auto" 
}

#for file in `ls 02_bt2_aligned/InPut*.bed | sort`; do
for file in `ls 02_bt2_aligned/H3K4me1.bed | sort`; do
        sample=`basename $file`
		sample=${sample:0:7}
        #submit1 $file $sample
		submit3.1 $file $sample
		#echo $file $sample
		#break
		#sleep 10
done

