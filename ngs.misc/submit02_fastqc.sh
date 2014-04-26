#!/bin/bash

submit_script() {
	qsub -b y -l h_vmem=8G -V -j y -m e -cwd -q all.q -N fastqc_$2 "fastqc -t \$NSLOTS -o 02_fastqc-sickle --noextract $1" 
}

#for folder in `ls -1 00_btangs/`; do
#	sample=${folder%????}
	#s1=00_btangs/$folder/${sample}_cleaned_a_a_1.fastq
	#s2=00_btangs/$folder/${sample}_cleaned_a_a_2.fastq
	#echo $sample $s1 $s2
	#break
	#submit_script $sample $s1
	# break
	#sleep 1
	#submit_script $sample $s2
	#sleep 1
	# break
#done

#for file in *.fastq.gz; do
#	echo ${file:0:6}
#done
#exit


#k=0
for file in `find 02_sickle -type f -name '*_merged.fq.gz'`; do
#	if [ $k = 0 ]
#	then
#		submit_script "Ep_1" $file
#		k+=1
#	else
#		submit_script "Ep_2" $file		
#	fi
	f1=`basename $file`
	f2=${f1%?????????????}
	# echo $f1 $f2
	submit_script $file $f2
	#sleep 1
done
