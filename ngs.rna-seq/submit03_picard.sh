#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -M mikhail.dozmorov@vcuhealth.org
#PBS -N picard
#PBS -j oe
# PBS -o /home/glasser/Dissertation/Dozmorov/RNA-seq/cuffdiffout3

cd $PBS_O_WORKDIR

# DIRIN=02_subread-align
DIRIN=03_sorted
DIROUT=04_picard

mkdir -p $DIROUT

for file in `find $DIRIN -type f -name '*.bam' | sort`; do
	sample=`basename $file .bam`;
	/home/mdozmorov/.local/bin/picard.jar CollectInsertSizeMetrics I=$file O=$DIROUT"/"$sample".txt" H=$DIROUT"/"$sample".pdf" M=0.5;
done

# picard.jar SortSam INPUT=$file OUTPUT=$DIROUT"/"$sample"_sorted-coord.bam" SORT_ORDER=coordinate
# picard.jar SortSam INPUT=$1 OUTPUT=$DIROUT"/"$2"_sorted-coord.bam" SORT_ORDER=coordinate 
# picard.jar MarkDuplicates INPUT=$DIROUT"/"$1"_sorted-coord.bam" OUTPUT=$DIROUT"/"$1"_sorted-coord_markdups.bam" METRICS_FILE=$DIROUT"/"$1"_sorted-coord_markdups.txt" REMOVE_DUPLICATES=true CREATE_INDEX=true
