#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -M mdozmorov@vcu.edu
#PBS -N samunmapped
#PBS -j oe
# PBS -o /home/glasser/Dissertation/Dozmorov/RNA-seq/cuffdiffout3

cd $PBS_O_WORKDIR

#FILE=/home/glasser/work.md/genes38_80.gtf

mkdir -p 03_unmapped

for file in `find 03_sorted/ -type f -name "*.bam" | sort`; do
        samtools view -f 4 $file > 03_unmapped/`basename $file .bam`.sam;
done
