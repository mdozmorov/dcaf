#!/bin/bash

submit_script() {
	qsub -b y -l h_vmem=8G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N sickle_$1 "sickle pe -f $2 -r $3 -t sanger -o 02_sickle/`basename $2` -p 02_sickle/`basename $3` -s 02_sickle/$1_single.fq" 
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


k=0
for file in `find 01_tmatic -type f -name '*_pair.fq.gz' | sort`; do
	if [ $k = 0 ]
	then
		f1=`basename $file`
		f2=$file
		k+=1
	else
		# echo ${f1%????????????????????????} $f1 $file
		submit_script ${f1%????????????????????????} $f2 $file	
		k=0
	fi
	# echo `basename $file`
	# submit_script $file `basename $file` # ${file:0:6}
	#sleep 10
done
