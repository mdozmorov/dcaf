#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=2
#PBS -M my@email.adr
#PBS -N cutadapt
#PBS -j oe
# PBS -o /path/to/stderr-stdout/output

cd $PBS_O_WORKDIR

# FASTA sequences of adapters
ADAPTERS=omicsoft.fa

# Input folder
DIRIN=00_raw/

# Output folder
DIROUT=01_trimmed

mkdir -p $DIROUT

# Paired end
k=0
for file in `find $DIRIN -name "*.fastq.gz" -type f | sort`; do
	if [ $k = 0 ]
	then
		s1=$file
		k+=1
	else
		s2=$file
		k=0
		cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o 01_cutadapt/`basename $s1` -p 01_cutadapt/`basename $s2` $s1 $s2;
	fi
done

# Single end
for file in `find $DIRIN -name "*.fastq.gz" -type f | grep "_sm_" | sort`; do
	# java -jar trimmomatic-0.33.jar SE -threads 12 -phred33 $file $DIROUT"/"`basename $file .fastq.gz`.fastq.gz ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36;
	cutadapt -a AAGATCGGAAGAGCACACGTCT -o $DIROUT"/"`basename $file` $file
done
