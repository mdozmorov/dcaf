#all:	.modules 00_raw-fqc 00_raw-fqc.txt 01_tmatic 01_tmatic-fqc 01_tmatic-fqc.txt 
#all:	02_star-raw
#all:	02_star-tmatic
#all:	02_bowtie2-raw
#all:	02_bowtie2-tmatic
#all:	.02_star-raw
#all:	.02_star-tmatic
#all:	.02_bowtie2-raw
#all:	.02_bowtie2-tmatic #Already included in 02_bowtie2-tmatic
#all:	02_star-raw-stats.txt
#all:	02_star-tmatic-stats.txt
#all:	02_bowtie2-raw-stats.txt
#all:	02_bowtie2-tmatic-stats.txt
#all:	03_markdup-star-raw
#all:	03_markdup-star-tmatic
#all:	03_markdup-bowtie2-raw
#all:	03_markdup-bowtie2-tmatic
#all:	.03_markdup-star-raw
#all:	.03_markdup-star-tmatic
#all:	.03_markdup-bowtie2-raw
#all:	.03_markdup-bowtie2-tmatic
all:	01_tmatic
#all:	03_markdup-star-raw-stats.txt
#all:	03_markdup-star-tmatic-stats.txt
#all:	03_markdup-bowtie2-raw-stats.txt
#all:	03_markdup-bowtie2-tmatic-stats.txt

SHELL=/bin/bash
#/usr/local/Modules/default/init/bash ; \

.modules:
	module add sge fastqc/0.9.4 sickle/1.21 cutadapt/1.2.1 trimmomatic/0.30 star/2.3.0e tophat/2.0.10 bowtie2/2.1.0 cufflinks/2.1.1 samtools/0.1.9 picard/1.99 htseq/0.5.3p9 macs/1.4.2 macs2/2.0.10.07132012 broadpeak/20130731 homer/4.5 sicer/1.1 ucsc/20110817 bedtools/2.19.1 python/3.4.0
	touch .modules

00_raw-fqc:
	LEN=22 ; \
	DIR="00_raw" ; \
	if [ ! -d $$DIR-fqc ]; then \
		mkdir $$DIR-fqc ; \
	fi ; \
	for file in `find $$DIR -type f -name '*.fq.gz'`; do \
		fname=`basename $$file` ; \
		fname=$${fname::$$LEN} ; \
		qsub -b y -l h_vmem=8G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-fqc "fastqc -o $$DIR-fqc --noextract $$file" ; \
	done

00_raw-fqc.txt:	00_raw-fqc
	python3 fastqc-summary 00_raw-fqc > 00_raw-fqc.txt

01_tmatic:	00_raw
	echo "tmatic exists"
	#DIR="01_tmatic" ; \
	#if [ ! -d $$DIR ]; then \
	#	mkdir $$DIR ; \
	#fi ; \
	#k=0 ; \
	#for file in `find 00_raw -type f -name *.fq.gz | sort`; do \
	#	if [ $$k = 0 ]; then \
	#		fname=`basename $$file` ; \
	#		f1=$$file ; \
	#		k+=1 ; \
	#	else \
	#		fname=$${fname::$$LEN} ; \
	#		f2=$$file ; \
	#		qsub -b y -pe threaded 4-24 -l h_vmem=8G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-tmatic "trimmomatic -threads \$$NSLOTS -phred33 $$f1 $$f2 $$DIR/`basename $$f1`_pair.fq.gz $$DIR/`basename $$f1`_unpair.fq.gz $$DIR/`basename $$f2`_pair.fq.gz $$DIR/`basename $$f2`_unpair.fq.gz CROP:70 ILLUMINACLIP:omicsoft.fa:4:40:10 LEADING:20 TRAILING:20 MINLEN:18" ; \
	#		k=0 ; \
	#	fi ; \
	#done

01_tmatic-fqc:	01_tmatic
	LEN=22 ; \
	DIR="01_tmatic" ; \
	if [ ! -d $$DIR-fqc ]; then \
		mkdir $$DIR-fqc ; \
	fi ; \
	for file in `find $$DIR -type f -name '*_pair.fq.gz'`; do \
		fname=`basename $$file` ; \
		fname=$${fname::$$LEN} ; \
		qsub -b y -l h_vmem=8G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-fqc "fastqc -o $$DIR-fqc --noextract $$file" ; \
	done

01_tmatic-fqc.txt:	01_tmatic-fqc
	python3 fastqc-summary 01_tmatic-fqc > 01_tmatic-fqc.txt

02_star-raw:	00_raw
	LEN=19 ; \
	DIR="02_star-raw" ; \
 	GDIR="/Volumes/hts_core/Shared/star_genomes/hg19/" ; \
	if [ ! -d $$DIR ]; then \
		mkdir $$DIR ; \
	fi ; \
	k=0 ; \
	for file in `find 00_raw -type f -name *.fq.gz | sort`; do \
		if [ $$k = 0 ]; then \
			fname=`basename $$file` ; \
			f1=$$file ; \
			k+=1 ; \
		else \
			fname=$${fname::$$LEN} ; \
			f2=$$file ; \
			k=0 ; \
			mkdir $$DIR/$$fname ; \
			qsub -b y -pe threaded 4-24 -l h_vmem=50G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-star-raw "STAR --genomeDir $$GDIR --readFilesIn $$f1 $$f2 --readFilesCommand zcat --runThreadN \$$NSLOTS --genomeLoad NoSharedMemory --outFileNamePrefix $$DIR/$$fname/" ; \
		fi ; \
	done

.02_star-raw:	02_star-raw
	DIR="02_star-raw" ; \
	LEN=19 ; \
	for file in `find $$DIR -type f -name '*.bam'`; do \
		fname=$${file:12:$$LEN} ; \
		qsub -b y -l h_vmem=20G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-star-raw2bam "picard SortSam INPUT=$$file OUTPUT=$$DIR/$$fname.csort.bam SORT_ORDER=coordinate CREATE_INDEX=true && rm $$file $$file.bai" ; \
		echo $$file $$fname ;\
	done
	#touch .02_star-raw
		#qsub -b y -l h_vmem=5G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-star-raw2bam "samtools view -bS $$file > $$DIR/$$fname.bam && picard SortSam INPUT=$$DIR/$$fname.bam OUTPUT=$$DIR/$$fname.coordsort.bam SORT_ORDER=coordinate CREATE_INDEX=true && rm $$DIR/$$fname.bam $$file" ; \

02_star-tmatic:	01_tmatic
	LEN=19 ; \
	DIR="02_star-tmatic" ; \
 	GDIR="/Volumes/hts_core/Shared/star_genomes/hg19/" ; \
	if [ ! -d $$DIR ]; then \
		mkdir $$DIR ; \
	fi ; \
	k=0 ; \
	for file in `find 01_tmatic -type f -name *_pair.fq.gz | sort`; do \
		if [ $$k = 0 ]; then \
			fname=`basename $$file` ; \
			f1=$$file ; \
			k+=1 ; \
		else \
			fname=$${fname::$$LEN} ; \
			f2=$$file ; \
			k=0 ; \
			mkdir $$DIR/$$fname ; \
			qsub -b y -pe threaded 4-24 -l h_vmem=50G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-star-tmatic "STAR --genomeDir $$GDIR --readFilesIn $$f1 $$f2 --readFilesCommand zcat --runThreadN \$$NSLOTS --genomeLoad NoSharedMemory --outFileNamePrefix $$DIR/$$fname/" ; \
		fi ; \
	done

.02_star-tmatic:	02_star-tmatic
	DIR="02_star-tmatic" ; \
	LEN=19 ; \
	for file in `find $$DIR -type f -name '*.bam'`; do \
		fname=$${file:15:$$LEN} ; \
		qsub -b y -l h_vmem=20G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-star-tmatic2bam "picard SortSam INPUT=$$file OUTPUT=$$DIR/$$fname.csort.bam SORT_ORDER=coordinate CREATE_INDEX=true && rm $$file $$file.bai" ; \
		echo $$file $$fname ;\
	done
	#touch .02_star-raw
		#qsub -b y -l h_vmem=5G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-star-raw2bam "samtools view -bS $$file > $$DIR/$$fname.bam && picard SortSam INPUT=$$DIR/$$fname.bam OUTPUT=$$DIR/$$fname.coordsort.bam SORT_ORDER=coordinate CREATE_INDEX=true && rm $$DIR/$$fname.bam $$file" ; \

02_bowtie2-raw:	00_raw
	LEN=19 ; \
	DIR="02_bowtie2-raw" ; \
 	GDIR="/Volumes/hts_core/Shared/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome" ; \
	if [ ! -d $$DIR ]; then \
		mkdir $$DIR ; \
	fi ; \
	k=0 ; \
	for file in `find 00_raw -type f -name *.fq.gz | sort`; do \
		if [ $$k = 0 ]; then \
			fname=`basename $$file` ; \
			f1=$$file ; \
			k+=1 ; \
		else \
			fname=$${fname::$$LEN} ; \
			f2=$$file ; \
			k=0 ; \
			qsub -b y -pe threaded 4-24 -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-bowtie2-raw "bowtie2 -p \$$NSLOTS --local --no-unal --mm -x $$GDIR -1 $$f1 -2 $$f2 | samtools view -bS - | samtools sort - $$DIR/$$fname.coordsort && samtools index $$DIR/$$fname.coordsort.bam" ; \
		fi ; \
	done

# Can include Picart SortSam into bowtie command
.02_bowtie2-raw:	02_bowtie2-raw
	DIR="02_bowtie2-raw" ; \
	LEN=19 ; \
	for file in `find $$DIR -type f -name '*.bam'`; do \
		fname=$${file:15:$$LEN} ; \
		echo $$file $$fname ;\
		qsub -b y -l h_vmem=20G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-bowtie2-raw2bam "picard SortSam INPUT=$$file OUTPUT=$$DIR/$$fname.csort.bam SORT_ORDER=coordinate CREATE_INDEX=true && rm $$file $$file.bai" ; \
	done
	touch .02_bowtie2-raw

02_bowtie2-tmatic:	01_tmatic
	LEN=19 ; \
	DIR="02_bowtie2-tmatic" ; \
 	GDIR="/Volumes/hts_core/Shared/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome" ; \
	if [ ! -d $$DIR ]; then \
		mkdir $$DIR ; \
	fi ; \
	k=0 ; \
	for file in `find 01_tmatic -type f -name *_pair.fq.gz | sort`; do \
		if [ $$k = 0 ]; then \
			fname=`basename $$file` ; \
			f1=$$file ; \
			k+=1 ; \
		else \
			fname=$${fname::$$LEN} ; \
			f2=$$file ; \
			k=0 ; \
			qsub -b y -pe threaded 4-24 -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-bowtie2-tmatic "bowtie2 -p \$$NSLOTS --local --no-unal --mm -x $$GDIR -1 $$f1 -2 $$f2 | samtools view -bS - > $$DIR/$$fname.bam && picard SortSam INPUT=$$DIR/$$fname.bam OUTPUT=$$DIR/$$fname.csort.bam SORT_ORDER=coordinate CREATE_INDEX=true && rm $$DIR/$$fname.bam" ; \
		fi ; \
	done

02_bowtie2-tmatic-stats.txt:	02_bowtie2-tmatic
	DIR="02_bowtie2-tmatic" ; \
	for file in $$DIR/*.bam; do \
		echo "Processing "$$file ; \
		samtools flagstat $$file >> $$DIR-stats.txt ; \
	done

03_markdup-star-raw:	02_star-raw
	LEN=19 ; \
	DIR="03_markdup-star-raw" ; \
	if [ ! -d $$DIR ]; then \
		mkdir $$DIR ; \
	fi ; \
	for file in `find 02_star-raw -type f -name '*.bam'`; do \
		fname=`basename $$file` ; \
		fname=$${fname::$$LEN} ; \
		qsub -b y -l h_vmem=8G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-markdup-star-raw "picard MarkDuplicates INPUT=$$file OUTPUT=$$DIR/$$fname.bam METRICS_FILE=$$DIR/$$fname.metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true" ; \
	done

# Can include Picart SortSam into bowtie command
.03_markdup-star-raw:	03_markdup-star-raw
	DIR="03_markdup-star-raw" ; \
	LEN=19 ; \
	for file in `find $$DIR -type f -name '*.bam'`; do \
		fname=$${file:20:$$LEN} ; \
		echo $$file $$fname ;\
		qsub -b y -l h_vmem=20G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-sort-markdup-star-raw "picard SortSam INPUT=$$file OUTPUT=$$DIR/$$fname.csort.bam SORT_ORDER=coordinate CREATE_INDEX=true && rm $$file $$file.bai" ; \
	done
	touch .03_markdup-star-raw

03_markdup-star-raw-stats.txt:	03_markdup-star-raw
	DIR="03_markdup-star-raw" ; \
	for file in $$DIR/*.bam; do \
		echo "Processing "$$file ; \
		samtools flagstat $$file >> $$DIR-stats.txt ; \
	done

03_markdup-star-tmatic:	02_star-tmatic
	LEN=19 ; \
	DIR="03_markdup-star-tmatic" ; \
	if [ ! -d $$DIR ]; then \
		mkdir $$DIR ; \
	fi ; \
	for file in `find 02_star-tmatic -type f -name '*.bam'`; do \
		fname=`basename $$file` ; \
		fname=$${fname::$$LEN} ; \
		qsub -b y -l h_vmem=8G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-markdup-star-tmatic "picard MarkDuplicates INPUT=$$file OUTPUT=$$DIR/$$fname.bam METRICS_FILE=$$DIR/$$fname.metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true" ; \
	done

# Can include Picart SortSam into bowtie command
.03_markdup-star-tmatic:	03_markdup-star-tmatic
	DIR="03_markdup-star-tmatic" ; \
	LEN=19 ; \
	for file in `find $$DIR -type f -name '*.bam'`; do \
		fname=$${file:23:$$LEN} ; \
		echo $$file $$fname ;\
		qsub -b y -l h_vmem=20G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-sort-markdup-star-tmatic "picard SortSam INPUT=$$file OUTPUT=$$DIR/$$fname.csort.bam SORT_ORDER=coordinate CREATE_INDEX=true && rm $$file $$file.bai" ; \
	done
	touch .03_markdup-star-tmatic

03_markdup-star-tmatic-stats.txt:	03_markdup-star-tmatic
	DIR="03_markdup-star-tmatic" ; \
	for file in $$DIR/*.bam; do \
		echo "Processing "$$file ; \
		samtools flagstat $$file >> $$DIR-stats.txt ; \
	done

03_markdup-bowtie2-raw:	02_bowtie2-raw
	LEN=19 ; \
	DIR="03_markdup-bowtie2-raw" ; \
	if [ ! -d $$DIR ]; then \
		mkdir $$DIR ; \
	fi ; \
	for file in `find 02_bowtie2-raw -type f -name '*.bam'`; do \
		fname=`basename $$file` ; \
		fname=$${fname::$$LEN} ; \
		qsub -b y -l h_vmem=8G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-markdup-bowtie2-raw "picard MarkDuplicates INPUT=$$file OUTPUT=$$DIR/$$fname.bam METRICS_FILE=$$DIR/$$fname.metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true" ; \
	done

# Can include Picart SortSam into bowtie command
.03_markdup-bowtie2-raw:	03_markdup-bowtie2-raw
	DIR="03_markdup-bowtie2-raw" ; \
	LEN=19 ; \
	for file in `find $$DIR -type f -name '*.bam'`; do \
		fname=$${file:23:$$LEN} ; \
		echo $$file $$fname ;\
		qsub -b y -l h_vmem=20G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-sort-markdup-bowtie2-raw "picard SortSam INPUT=$$file OUTPUT=$$DIR/$$fname.csort.bam SORT_ORDER=coordinate CREATE_INDEX=true && rm $$file $$file.bai" ; \
	done
	touch .03_markdup-bowtie2-raw

03_markdup-bowtie2-raw-stats.txt:	03_markdup-bowtie2-raw
	DIR="03_markdup-bowtie2-raw" ; \
	for file in $$DIR/*.bam; do \
		echo "Processing "$$file ; \
		samtools flagstat $$file >> $$DIR-stats.txt ; \
	done

03_markdup-bowtie2-tmatic:	02_bowtie2-tmatic
	LEN=19 ; \
	DIR="03_markdup-bowtie2-tmatic" ; \
	if [ ! -d $$DIR ]; then \
		mkdir $$DIR ; \
	fi ; \
	for file in `find 02_bowtie2-tmatic -type f -name '*.bam'`; do \
		fname=`basename $$file` ; \
		fname=$${fname::$$LEN} ; \
		qsub -b y -l h_vmem=8G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-markdup-bowtie2-tmatic "picard MarkDuplicates INPUT=$$file OUTPUT=$$DIR/$$fname.bam METRICS_FILE=$$DIR/$$fname.metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true" ; \
	done

# Can include Picart SortSam into bowtie command
.03_markdup-bowtie2-tmatic:	03_markdup-bowtie2-tmatic
	DIR="03_markdup-bowtie2-tmatic" ; \
	LEN=19 ; \
	for file in `find $$DIR -type f -name '*.bam'`; do \
		fname=$${file:26:$$LEN} ; \
		echo $$file $$fname ;\
		qsub -b y -l h_vmem=20G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N job-$$fname-sort-markdup-bowtie2-tmatic "picard SortSam INPUT=$$file OUTPUT=$$DIR/$$fname.csort.bam SORT_ORDER=coordinate CREATE_INDEX=true && rm $$file $$file.bai" ; \
	done
	touch .03_markdup-bowtie2-tmatic

03_markdup-bowtie2-tmatic-stats.txt:	03_markdup-bowtie2-tmatic
	DIR="03_markdup-bowtie2-tmatic" ; \
	for file in $$DIR/*.bam; do \
		echo "Processing "$$file ; \
		samtools flagstat $$file >> $$DIR-stats.txt ; \
	done



clean:
	-rm *fqc.o*
	-rm *.po*
	-mv job*tmatic.o* jobs/
