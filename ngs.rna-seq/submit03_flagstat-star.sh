#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -M my@email.adr
#PBS -N flagstat
#PBS -j oe
# PBS -o /path/to/stderr-stdout/output

cd $PBS_O_WORKDIR

# Input folder
DIRIN=02_star-align/

# Output file
FOUT=flagstat_star-align.txt

for file in `find $DIRIN -type f -name "*.bam" | sort`; do
		echo $file >> $FOUT;
        samtools flagstat $file >> $FOUT;
done
