#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=10
#PBS -M mdozmorov@vcu.edu
#PBS -N cuffquant
#PBS -j oe
# PBS -o /home/glasser/Dissertation/Dozmorov/RNA-seq/cuffdiffout3

cd $PBS_O_WORKDIR

cufflinks -p 10 -o 04_cufflinks -G /home/sequencing/data/ExtData/UCSC/hg38/Homo_sapiens.GRCh38.83.gtf 03_sorted/175278A_S6.bam 03_sorted/182696A_S8.bam 
# 03_sorted/181040A_S5.bam 03_sorted/NSG21A_S7.bam 03_sorted/180092B_S3.bam 03_sorted/183291B_S4.bam 03_sorted/178200B_S1.bam 03_sorted/180090B_S2.bam 03_sorted/182835C_S10.bam 03_sorted/183303C_S11.bam 03_sorted/181658C_S9.bam 03_sorted/182837C_S12.bam

for file in 03_sorted/*.bam; do 
	echo $file;
	cat > submit04_cuffnorm_`basename $file .bam`.sh << EOT
#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=10
#PBS -M mdozmorov@vcu.edu
#PBS -N cufflinks
#PBS -j oe
# PBS -o /home/glasser/Dissertation/Dozmorov/RNA-seq/cuffdiffout3
cd \$PBS_O_WORKDIR
cufflinks -p 10 -o 04_cufflinks/`basename $file .bam` -G /home/sequencing/data/ExtData/UCSC/hg38/Homo_sapiens.GRCh38.83.gtf $file
EOT

qsub submit04_cuffnorm_`basename $file .bam`.sh;

done
