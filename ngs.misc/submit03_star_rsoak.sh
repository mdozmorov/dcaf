#!/bin/bash

submit_script() {
	qsub -b y -pe threaded 4-24 -l h_vmem=50G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N star_rsoak_$1 "mkdir 03_star_rsoak/$1 && STAR --genomeDir ../../RepeatSoaker/index.star --readFilesIn $2 $3 --runThreadN \$NSLOTS --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 --outFilterMatchNmin 10 --genomeLoad NoSharedMemory --outFileNamePrefix 03_star_rsoak/$1/ --outReadsUnmapped Fastx" 
}
# --genomeLoad LoadAndRemove
k=0
for file in `find 02_sickle -type f -name *_pair.fq.gz | sort`; do
	if [ $k = 0 ]
	then
		f1=`basename $file`
		f2=$file
		k+=1
	else
		echo ${f1%????????????????????} $f1 $file
		submit_script ${f1%????????????????????} $f2 $file	
		k=0
	fi
	# echo `basename $file`
	# submit_script $file `basename $file` # ${file:0:6}
	sleep 3
done