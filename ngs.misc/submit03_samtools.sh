#!/bin/bash

submit1bt() {
	qsub -b y -pe threaded 8 -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $3_bt2 "bowtie2 -p 8 --local --very-sensitive-local --no-unal --mm -x /Volumes/hts_core/Shared/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome -1 $1 -2 $2 -S aligned_cleaned/$3" 
}


submit1() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $2_st "samtools view -bS $1 | samtools sort - aligned_sorted/$2.sorted && samtools index aligned_sorted/$2.sorted.bam" 
}

submit2() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $2_s2b "samtools view -bS $1 > $2.bam && rm $1" 
}

submit3() {
	qsub -b y -l h_vmem=5G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $2_sidx "samtools index $1" 
}

for file in `ls aligned_raw/*.bam | sort`; do
                sample=`basename $file`
		sample=${sample%????????}
                submit3 $file $sample
		echo $file $sample
		#break
		sleep 10
done

