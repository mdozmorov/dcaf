#!/bin/bash
submit1() {
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N multicov_macs2 "bedtools multicov -bams $k -bed broad_peak_coords.s.m.bed > broad_peak_mtx.txt"
}
k=""
for file in `find 03_bowtie-tmatic -type f -name '*.bam' | sort`; do
	k=$k" "$file
	#echo $file
done
bedtools multicov -bams $k -bed macs2-bowtie_peaks.m.s.bed  > macs2-bowtie_peaks_mtx.txt
