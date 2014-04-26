#!/bin/bash

submit1bt() {
	qsub -b y -pe threaded 8 -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $3_bt2 "bowtie2 -p 8 --local --very-sensitive-local --no-unal --mm -x /Volumes/hts_core/Shared/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome -1 $1 -2 $2 -S aligned_cleaned/$3" 
}


submit1() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $2_pfS "picard SortSam INPUT=$1 OUTPUT=02_aligned/Controls/$2.s.bam SORT_ORDER=coordinate && rm $1" 
}

submit2() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $2_pMD "picard MarkDuplicates INPUT=$1 OUTPUT=03_aligned_picard/Controls/$2.f.bam METRICS_FILE=03_aligned_picard/Controls/$2.metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true" 
}

submit12() {
	qsub -b y -pe threaded 10-24 -l h_vmem=50G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $3_bt1 "bowtie -p \$NSLOTS --chunkmbs 512 -S --trim3 5 /Volumes/hts_core/Shared/bowtie_index/mm9/mm9 -1 $1 -2 $2 | samtools view -bS - | samtools sort - aligned2/$3.sorted && samtools index aligned2/$3.sorted.bam" 
}

submit3() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $2_pHead "picard AddOrReplaceReadGroups INPUT=$1 OUTPUT=$3/$2.h.bam SORT_ORDER=coordinate RGID=Chris RGLB=MS RGPL=Illumina RGSM=$2 RGPU=$2 CREATE_INDEX=true && rm $1" 
}

for file in `find 03_aligned_picard/ -type f -name '*.bam'`; do
                sample=`basename $file`
				basedir=$(dirname $file)
		sample=${sample:0:${#sample}-4}
                submit3 $file $sample $basedir
		#echo $file $sample $basedir
		#break
		sleep 1
done

