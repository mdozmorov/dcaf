#!/bin/bash

submit1() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $2_htseq_ensembl "samtools view -h $1 | htseq-count -q - /home/dozmorovm/work/Homo_sapiens.GRCh37.69.cleaned.gtf > aligned_raw/$2.counts.ensembl.txt"
}

for file in `ls aligned_raw/*.bam | sort`; do
	sample=`basename $file`
	sample=${sample%????????}
	# sd=${ss}accepted.bam
	echo $file $sample
	submit1 $file $sample
	# break
	sleep 10
done
