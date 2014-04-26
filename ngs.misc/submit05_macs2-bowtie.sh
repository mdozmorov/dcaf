#!/bin/bash

submit1() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N macs2_$2 "macs2 callpeak -t $1 -B --name 05_macs2-bowtie/$2/$2" 
}

submit2() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N macs2_$2 "macs2 callpeak -t $1 -c 02_bt2_aligned/InPut_D.sorted.bam -B --broad --name $2 --shift-control --nolambda" 
}

submit3() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N s03_m2_$2 "macs2 callpeak -t $1 -B --broad --name $2" 
}

for file in `find 03_bowtie-tmatic -type f -name '*.bam' | sort`; do
        sample=`basename $file`
		#sample=`echo $sample | tr "/" "\t" | cut -f2`
        mkdir 05_macs2-bowtie/${sample%??????????????}
		#echo $file ${sample%??????????????}		
		submit1 $file ${sample%??????????????}	
done

