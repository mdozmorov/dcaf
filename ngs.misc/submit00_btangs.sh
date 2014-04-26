#!/bin/bash

BOWTIE_IDX=/Volumes/hts_core/Shared/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome


submit_script() {
	qsub -b y -l h_vmem=8G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $1_00_btangs "clean_sample.rb -s $1 -r a -l a -b 00_btangs/ $2 $3" 
}

SAMPLE_DIR=/Volumes/hts_raw/runs/120518_SN470_0220_BC0VMEACXX/Analysis/rnaseq_samples/unaligned/Project_Chris_RNA_SEQ/
#SAMPLE_DIR=/Volumes/hts_raw/runs/120517_SN470_0219_AD12NRACXX/Analysis/rnaseq_samples/unaligned/Project_Chris_RNA_SEQ/

# submit_script test $1 $2
# echo ${1##*/}
# exit

for folder in $SAMPLE_DIR/*; do
	sample=`basename $folder`
	sample=${sample:7:100}
	s1=`ls $folder/*.fastq.gz | grep -P 'L00\d_R1'`
	s2=`ls $folder/*.fastq.gz | grep -P 'L00\d_R2'`
	echo $sample
	echo $s1
	echo $s2
		
	submit_script $sample $s1 $s2
	# break
	sleep 1
	# submit_script $sample $s2
	# sleep 1
	# break
done
