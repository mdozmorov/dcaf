#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=24
#PBS -M my@email.adr
#PBS -N starindex
#PBS -j oe
# PBS -o /path/to/stderr-stdout/output

# Path to genome annotation files
# Note PBS jobs have no internet access. Download genome annotation files manually.
# Make $DIRINDEX folder
# Download genome annotation files from https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files
# wget https://gdc-api.nci.nih.gov/data/62f23fad-0f24-43fb-8844-990d531947cf; mv 62f23fad-0f24-43fb-8844-990d531947cf GRCh38.d1.vd1.fa.tar.gz; tar -zxvf GRCh38.d1.vd1.fa.tar.gz
# wget https://gdc-api.nci.nih.gov/data/fe1750e4-fc2d-4a2c-ba21-5fc969a24f27; mv fe1750e4-fc2d-4a2c-ba21-5fc969a24f27 gencode.v22.annotation.gtf.gz; gzip -d gencode.v22.annotation.gtf.gz
DIRINDEX=/home/mdozmorov/sequencing/data/ExtData/UCSC/hg38gdc
INDEXGENOME=$DIRINDEX"/GRCh38.d1.vd1.fa"
INDEXGENE=$DIRINDEX"/gencode.v22.annotation.gtf"

# DIRINDEX=/home/sequencing/data/ExtData/UCSC/mm10
# INDEXGENOME=$DIRINDEX"/mm10.fa"
# INDEXGENE=$DIRINDEX"/Mus_musculus.GRCm38.83_filt.gtf"


cd $PBS_O_WORKDIR

STAR --runThreadN 23 --runMode genomeGenerate --genomeDir $DIRINDEX --genomeFastaFiles $INDEXGENOME --sjdbGTFfile $INDEXGENE
