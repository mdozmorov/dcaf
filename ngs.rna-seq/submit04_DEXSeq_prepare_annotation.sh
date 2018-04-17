#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=2
#PBS -M my@email.adr
#PBS -N DEXSeq_annot
#PBS -j oe
# PBS -o /path/to/stderr-stdout/output

cd $PBS_O_WORKDIR

# DEXSeq path to Python scripts
DIRDEXSeq=/home/mdozmorov/R/x86_64-redhat-linux-gnu-library/3.4/DEXSeq/python_scripts
# dexseq_count.py  dexseq_prepare_annotation.py

# Human settings
# DIRGENOME=/home/sequencing/data/ExtData/UCSC/hg38gdc
# GFTIN=gencode.v22.annotation.gtf
# GTFOUT=gencode.v22.annotation_DEXSeq.gff

# Mouse settings
DIRGENOME=/home/sequencing/data/ExtData/UCSC/mm10
GFTIN=Mus_musculus.GRCm38.83_filt.gtf
GTFOUT=Mus_musculus.GRCm38.83_filt_DEXSeq.gff

python $DIRDEXSeq"/dexseq_prepare_annotation.py" $DIRGENOME"/"$GFTIN $DIRGENOME"/"$GTFOUT