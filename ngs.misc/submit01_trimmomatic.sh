#!/bin/bash

submit_script() {
	qsub -b y -pe threaded 4-24 -l h_vmem=8G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N tmatic_$1 "trimmomatic -threads \$NSLOTS -phred33 $2 $3 $2_pair.fq.gz $2_unpair.fq.gz $3_pair.fq.gz $3_unpair.fq.gz CROP:70 ILLUMINACLIP:omicsoft.fa:4:40:10 LEADING:20 TRAILING:20 MINLEN:18" 
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
for file in `find raw -type f -name *.fq.gz | sort`; do
	if [ $k = 0 ]
	then
		f1=`basename $file`
		f2=$file
		k+=1
	else
		#echo ${f1:3:6} $f1 $file
		#echo ${f1%?????????????} $f2 $file	
		submit_script ${f1%?????????????} $f2 $file	
		k=0
	fi
	# echo `basename $file`
	# submit_script $file `basename $file` # ${file:0:6}
	#sleep 10
done
