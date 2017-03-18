#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=4
#PBS -M my@email.adr
#PBS -N samsort
#PBS -j oe
# PBS -o /path/to/stderr-stdout/output

cd $PBS_O_WORKDIR

# Input folder
DIRIN=02_subread-align/

# Output folder
DIROUT=03_sorted

mkdir -p $DIROUT

for file in `find $DIRIN -type f -name "*.bam" | sort`; do
        samtools sort $file $DIROUT"/"`basename $file .bam`;
done
