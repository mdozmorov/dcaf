#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=24
#PBS -l mem=24000MB
#PBS -M my@email.adr
#PBS -N featureCounts
#PBS -j oe
# PBS -o /path/to/stderr-stdout/output

cd $PBS_O_WORKDIR

# Path to genome annotation files
DIRINDEX=/path/to/genome/annotation/files
INDEXGENE=gencode.v22.annotation.gtf.gz

# Output file with gene counts
COUNTS=counts.txt

# Rscript rscript_featureCount.R
# featureCounts -p -t exon -g gene_id -a /home/sequencing/data/ExtData/UCSC/hg38/Homo_sapiens.GRCh38.83.gtf -o counts.txt 02_subread-align/AI1_ATCACG_R1.bam 02_subread-align/AI2_ACTTGA_R1.bam 02_subread-align/AI3_GTTTCG_R1.bam 02_subread-align/AI4_ACAGTG_R1.bam 02_subread-align/AI5.bam 02_subread-align/AI6.bam 02_subread-align/AI7_ACTGAT_R1.bam 02_subread-align/AI8_TGACCA_R1.bam 02_subread-align/AI9_TAGCTT_R1.bam 02_subread-align/AI10_ATGTCA_R1.bam 02_subread-align/AI11_CGATGT_R1.bam 02_subread-align/AI12_TTAGGC_R1.bam 02_subread-align/AI13_GATCAG_R1.bam 02_subread-align/AI14_AGTCAA_R1.bam 

# Summarize counts per gene_id.
# For single-end, remove "-p" option
featureCounts -T 24 -p -t exon -g gene_id -a $DIRINDEX"/"$INDEXGENE -o $COUNTS <comma-separated relative paths to BAM files>
