#!/bin/bash

submit_script() {
	qsub -b y -l h_vmem=8G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N fastqc_$2 "fastqc -t \$NSLOTS -o 00_fastqc --noextract $1" 
}

for folder in `find raw -type d`; do
	if [ $folder = "raw" ] # Skip the first folder returned by find, as it is just root folder
	then
		continue
	fi
	f=`find $folder -type f -name '*.fastq.gz' | sort` # Get the lsit of files in the current folder
	flist=( $f ) # Convert to array
	fname=`basename ${flist[0]}` # Get the file name
	fname=${fname%????????????????????} # Trim it
	zcat ${flist[0]} ${flist[2]} | gzip > "raw/"$fname"R1_merged.fq.gz"
	zcat ${flist[1]} ${flist[3]} | gzip > "raw/"$fname"R2_merged.fq.gz"
done
