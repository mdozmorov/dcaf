#!/bin/bash

if [ ! -f hg19.genome ]; then
	mysql -ugenome -hgenome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo" > hg19.genome
fi

function shuffle () {
	base="chr1\t1\t$((1+$1))"
	for i in `seq 1 $2`; do
		echo -e $base >> $$tmp.bed
	done
	shuffleBed -g hg19.genome -i $$tmp.bed
	rm -f $$tmp.bed
}

N_BP=( 1 10 100 1000 )
N_LINES=( 10 100 1000 )

for bp in "${N_BP[@]}"; do
	mkdir -p random/$bp
	for lines in "${N_LINES[@]}"; do
		shuffle $bp $lines > random/$bp/$lines.bed
	done
done
#rm $$tmp.bed


