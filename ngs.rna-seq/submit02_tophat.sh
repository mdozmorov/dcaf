mkdir -p 02_tophat_mm9

# # Paired end
# k=0
# for file in `find 01_trimmed/ -name "*_paired.fastq.gz" -type f | sort`; do
# 	if [ $k = 0 ]
# 	then
# 		s1=$file
# 		k+=1
# 	else
# 		s2=$file
# 		k=0
# 		cat > submit02_tophat_`basename $s1 _L001_R1_001_paired.fastq.gz`.sh << EOT
# #!/bin/bash
# #PBS -S /bin/bash
# #PBS -V
# #PBS -l nodes=1:ppn=10
# #PBS -M mdozmorov@vcu.edu
# #PBS -N tophat
# #PBS -j oe
# # PBS -o /home/glasser/Dissertation/Dozmorov/RNA-seq/cuffdiffout3
# cd \$PBS_O_WORKDIR
# tophat -p 10 -o 02_tophat/`basename $s1 _L001_R1_001_paired.fastq.gz` -G /home/sequencing/data/ExtData/UCSC/mm10/Mus_musculus.GRCm38.83_filt.gtf /home/sequencing/data/ExtData/UCSC/mm10/mm10 $s1 $s2
# EOT

# qsub submit02_tophat_`basename $s1 _L001_R1_001_paired.fastq.gz`.sh;
# 	fi
# done

# Single end
for file in `find 01_trimmed/ -name "*.fastq.gz" -type f | sort | tail -2`; do
	cat > submit02_tophat_mm9_`basename $file _L002_R1_001.fastq.gz`.sh << EOT
#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=2
#PBS -M mdozmorov@vcu.edu
#PBS -N tophat
#PBS -j oe
# PBS -o /home/glasser/Dissertation/Dozmorov/RNA-seq/cuffdiffout3
SHELL=/bin/tcsh
cd \$PBS_O_WORKDIR
# bowtie2 -p 10 -x /home/sequencing/data/ExtData/UCSC/mm10/mm10 -U $file | samtools view -bS - | samtools sort - 02_bowtie2/`basename $file .fastq.gz` && samtools index 02_bowtie2/`basename $file .fastq.gz`.bam
tophat -o 02_tophat_mm9/`basename $file _L002_R1_001.fastq.gz` -G /home/sequencing/data/ExtData/UCSC/mm9/Mus_musculus.NCBIM37.67_filt.gtf /home/sequencing/data/ExtData/UCSC/mm9/mm9 $file
EOT

qsub submit02_tophat_mm9_`basename $file _L002_R1_001.fastq.gz`.sh;

done

# bowtie2 -p 10 -x /home/sequencing/data/ExtData/UCSC/mm10/mm10 -U $file | samtools view -bS - | samtools sort - 02_bowtie2/`basename $file .fastq.gz` && samtools index 02_bowtie2/`basename $file .fastq.gz`.bam
