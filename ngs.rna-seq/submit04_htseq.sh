#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=2
#PBS -M mdozmorov@vcu.edu
#PBS -N htseq
#PBS -j oe
# PBS -o /home/glasser/Dissertation/Dozmorov/RNA-seq/cuffdiffout3

cd $PBS_O_WORKDIR

FILEGTF=/home/sequencing/data/ExtData/UCSC/rn6/Rattus_norvegicus.Rnor_6.0.85_chr.gtf

mkdir 04_htseq

for file in 03_sorted/*; do
        echo $file;
        k=`basename $file .bam`;
        htseq-count -f bam -r name -s no $file $FILEGTF > 04_htseq/$k.txt;
done
