#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=24
#PBS -M mdozmorov@vcu.edu
#PBS -N samsortstar
#PBS -j oe
# PBS -o /home/glasser/Dissertation/Dozmorov/RNA-seq/cuffdiffout3

cd $PBS_O_WORKDIR

#FILE=/home/glasser/work.md/genes38_80.gtf

# Sort by coordinate
mkdir -p 03_sorted-star

for file in `find 02_star-align/ -type f -name "*.bam" | sort`; do
	/home/mdozmorov/.local/bin/samtools sort -@ 23 -o 03_sorted-star/`basename $file` -O bam $file;
done

for file in `find 03_sorted-star/ -type f -name "*.bam" | sort`; do
	/home/mdozmorov/.local/bin/samtools index -@ 23 $file;
done

# Sort by name
mkdir -p 03_sorted-star.n

for file in `find 02_star-align/ -type f -name "*.bam" | sort`; do
	/home/mdozmorov/.local/bin/samtools sort -@ 23 -n -o 03_sorted-star.n/`basename $file` -O bam $file;
done

