#!/bin/bash

submit1bt() {
	qsub -b y -pe threaded 8 -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $3_bt2 "bowtie2 -p 8 --local --very-sensitive-local --no-unal --mm -x /Volumes/hts_core/Shared/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome -1 $1 -2 $2 -S aligned_cleaned/$3" 
}


submit1() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N s02_st_$2 "samtools view -bS $1 | samtools sort - 02_bt2_aligned/$2.sorted && samtools index 02_bt2_aligned/$2.sorted.bam && rm $1" 
}

submit2() {
	qsub -b y -l h_vmem=5G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N samtools_star_$2 "samtools view -bS $1 > 03_star_bam/$2.bam && samtools index 03_star_bam/$02.bam" 
}

submit3() {
	qsub -b y -l h_vmem=5G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $2_sidx "samtools index $1" 
}

for file in `ls 02_bt2_aligned/*.sam | sort`; do
        sample=`basename $file`
		sample=${sample:0:7}
        submit1 $file $sample
		#echo $file $sample
		#break
		#sleep 10
done

for file in `find 03_star -type f -name '*.sam' | sort`; do
#for file in `ls 02_bt2_aligned/PU*.bam | sort`; do
        sample=`dirname $file`
		sample=`echo $sample | tr "/" "\t" | cut -f2`
        #mkdir 05_macs2/$sample
		submit2 $file $sample
		# submit2 $file $sample
		#echo $file $sample
		#break
		#sleep 10
done
