#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=12
#PBS -M my@email.adr
#PBS -N subread
#PBS -j oe
# PBS -o /path/to/stderr-stdout/output

cd $PBS_O_WORKDIR

# Input folder
DIRIN=01_trimmed/
# Output folder
DIROUT=02_subread-align
# Path to genome annotation files
DIRINDEX=/path/to/genome/annotation/files
# DIRINDEX=/home/mdozmorov/sequencing/data/ExtData/UCSC/hg38gdc

mkdir -p $DIROUT

# Paired end
# k=0
# for file in `find $DIRIN -name "*_paired.fastq.gz" -type f | sort`; do
# 	if [ $k = 0 ]
# 	then
# 		s1=$file
# 		k+=1
# 	else
# 		s2=$file
# 		k=0
# 		subread-align -T 12 --gzFASTQinput --BAMoutput -i $DIRINDEX -r $s1 -R $s2 -o $DIROUT"/"`basename $s1 _L001_R1_001_paired.fastq.gz`.bam;
# 	fi
# done

for file in `find $DIRIN -name "*.fastq.gz" -type f | sort`; do
	subread-align -T 12 --gzFASTQinput --BAMoutput -i $DIRINDEX -r $file -o $DIROUT"/"`basename $file .fastq.gz`.bam;
done
