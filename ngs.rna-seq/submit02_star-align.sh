#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=24
#PBS -M my@email.adr
#PBS -N star
#PBS -j oe
# PBS -o /path/to/stderr-stdout/output

cd $PBS_O_WORKDIR

# Input folder
DIRIN=01_trimmed/
# Output folder
DIROUT=02_star-align
# Path to genome annotation files
# DIRINDEX=/home/sequencing/data/ExtData/UCSC/mm10
DIRINDEX=/home/mdozmorov/sequencing/data/ExtData/UCSC/hg38gdc

mkdir -p $DIROUT

# Paired end
k=0
for file in `find $DIRIN -name "*_paired.fastq.gz" -type f | grep -v Undetermined | sort`; do
	if [ $k = 0 ]
	then
		s1=$file
		k+=1
		prefix=`basename $file _R1_001_paired.fastq.gz`;
	else
		s2=$file
		k=0
		STAR  --runThreadN 23 --genomeDir $DIRINDEX --readFilesIn $s1 $s2 --readFilesCommand zcat --outSAMtype BAM Unsorted --outSAMorder Paired --outReadsUnmapped Fastx --outFileNamePrefix $DIROUT/$prefix. --outFilterMultimapNmax 10;
	fi
done

# Single end
# for file in `find $DIRIN -name "*.fastq.gz" -type f | sort`; do
# 	prefix=`basename $file .fastq.gz`;
# 	STAR  --runThreadN 23 --genomeDir $DIRINDEX --readFilesIn $file --readFilesCommand zcat --outSAMtype BAM Unsorted --outSAMorder Paired --outReadsUnmapped Fastx --outFileNamePrefix $DIROUT/$prefix. --outFilterMultimapNmax 10;
# done
