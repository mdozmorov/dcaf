#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=24
#PBS -M mdozmorov@vcu.edu
#PBS -N mdozsamsort
#PBS -j oe
# PBS -o /home/glasser/Dissertation/Dozmorov/RNA-seq/cuffdiffout3

cd $PBS_O_WORKDIR

#FILE=/home/glasser/work.md/genes38_80.gtf

mkdir 03_sorted


for file in `find 02_subread-align/ -type f -name "*.bam" | sort`; do
	/home/mdozmorov/.local/bin/samtools sort -@ 23 -o 03_sorted/$file -O BAM $file;
done

for file in `find 03_sorted/ -type f -name "*.bam" | sort`; do
	/home/mdozmorov/.local/bin/samtools index $file;
done
