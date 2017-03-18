#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=12
#PBS -M my@email.adr
#PBS -N fastqcbam
#PBS -j oe
# PBS -o /path/to/stderr-stdout/output

cd $PBS_O_WORKDIR

# Input folder
DIRIN=03_sorted/

# Output folder
DIROUT=03_fastqc-bam

mkdir -p $DIROUT

for file in `find $DIRIN -type f -name "*.bam" | sort`; do 
	fastqc -t 12 -o $DIROUT --noextract $file; 
done

# Extract tab-separated FASTQC summary
python fastqc-summary -s $DIROUT > $DIROUT".txt"