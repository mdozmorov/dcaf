#!/bin/bash

submit1() {
	qsub -b y -pe threaded 4-12 -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $2_cl "cufflinks -p \$NSLOTS -g /home/dozmorovm/work/hg19.knowngene.gtf -M /home/dozmorovm/work/hg19.knowngene.mask.gtf -b /Volumes/hts_core/Shared/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome.fa -u --library-type=fr-secondstrand -o cl_out/$2 $1"
}

for file in `ls aligned_sorted/*.bam`; do
	sample=`basename $file`
	sample=${sample%???????????}
	# sd=${ss}accepted.bam
	echo $file $sample
	submit1 $file $sample
	# break
	sleep 10
done
