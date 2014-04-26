#!/bin/bash

submit1bt() {
	qsub -b y -pe threaded 5-12 -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N bt2_$3 "bowtie2 -p \$NSLOTS  --local --very-sensitive-local --no-unal --mm -x /Volumes/hts_core/Shared/bowtie2_index/mm9/mm9 -1 $1 -2 $2 | samtools view -bS - > $3.bam" 
}


submit1() {
	qsub -b y -pe threaded 10-24 -l h_vmem=50G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N bt2_$3 "bowtie2 -p \$NSLOTS --local --very-fast-local --no-unal --mm -x /Volumes/hts_core/Shared/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome -1 $1 -2 $2 | samtools view -bS - | samtools sort - 03_bowtie-tmatic/$3.sorted && samtools index 03_bowtie-tmatic/$3.sorted.bam" 
}

ubmit11() {
	qsub -b y -pe threaded 10-24 -l h_vmem=50G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $3_bt2 "bowtie2 --trim3 5 --phred33 --very-sensitive --no-unal --mm -x /Volumes/hts_core/Shared/bowtie2_index/mm9/mm9 -1 $1 -2 $2 | samtools view -bS - | samtools sort - aligned1/$3.sorted && samtools index aligned1/$3.sorted.bam" 
}

submit12() {
	qsub -b y -pe threaded 10-24 -l h_vmem=50G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $3_bt1 "bowtie -p \$NSLOTS --chunkmbs 512 -S --trim3 5 /Volumes/hts_core/Shared/bowtie_index/mm9/mm9 -1 $1 -2 $2 | samtools view -bS - | samtools sort - aligned2/$3.sorted && samtools index aligned2/$3.sorted.bam" 
}

k=0
for file in `find 01_tmatic -type f -name *_pair.fq.gz | sort`; do
	if [ $k = 0 ]
	then
		f1=`basename $file`
		f2=$file
		k+=1
	else
		#echo ${f1%????????????????????????} $f1 $file
		submit1 $f2 $file ${f1%????????????????????????} 
		k=0
	fi
done

