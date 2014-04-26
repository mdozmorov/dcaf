#!/bin/bash

submit1() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N $2_htseq_ensembl "samtools pileup $1 | perl -ne `'BEGIN{print "track type=wiggle_0 name=$2 description=$2\n"};($c, $start, undef, $depth) = split; if ($c ne $lastC) { print "variableStep chrom=$c\n"; };$lastC=$c;next unless $. % 10 ==0;print "$start\t$depth\n" unless $depth<3;' > $1%????.wig`"
}

for file in `find . -name *rd.bam | sort`; do
	sample=`basename $file`
	echo $file $sample
	# submit1 $file $sample
	# break	
done
