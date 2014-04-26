#all:	04_rs-markdup-bowtie2-tmatic


04_rs-bowtie2-raw:	02_bowtie2-raw
	LEN=19 ; \
	GDIR="/Volumes/hts_core/Shared/repeat_soaker/rmsk.hg19.bed" ; \
	DIR="04_rs-bowtie2-raw" ; \
	limits="0.75 0.50 0.25 0.01 0.00" ; \
	if [ ! -d $$DIR ]; then \
		mkdir $$DIR ; \
	fi ; \
	for file in `find 02_bowtie2-raw -type f -name '*.bam'`; do \
		fname=`basename $$file` ; \
		fname=$${fname%??????????} ; \
		for i in $$limits; do \
			qsub -b y -l h_vmem=8G -V -j y -m e -cwd -q all.q -N job-$$fname-rs$$i-bowtie2-raw "repeat-soaker -r $$GDIR -p $$i -o $$DIR/$$fname.$$i.bam $$file" ; \
		done ; \
    done

04_rs-bowtie2-tmatic:	02_bowtie2-tmatic
	LEN=19 ; \
	GDIR="/Volumes/hts_core/Shared/repeat_soaker/rmsk.hg19.bed" ; \
	DIR="04_rs-bowtie2-tmatic" ; \
	limits="0.75 0.50 0.25 0.01 0.00" ; \
	if [ ! -d $$DIR ]; then \
		mkdir $$DIR ; \
	fi ; \
	for file in `find 02_bowtie2-tmatic -type f -name '*.bam'`; do \
		fname=`basename $$file` ; \
		fname=$${fname%??????????} ; \
		for i in $$limits; do \
			qsub -b y -l h_vmem=8G -V -j y -m e -cwd -q all.q -N job-$$fname-rs$$i-bowtie2-tmatic "repeat-soaker -r $$GDIR -p $$i -o $$DIR/$$fname.$$i.bam $$file" ; \
		done ; \
    done


04_rs-markdup-bowtie2-raw:	03_markdup-bowtie2-raw
	LEN=19 ; \
	GDIR="/Volumes/hts_core/Shared/repeat_soaker/rmsk.hg19.bed" ; \
	DIR="04_rs-markdup-bowtie2-raw" ; \
	limits="0.75 0.50 0.25 0.01 0.00" ; \
	if [ ! -d $$DIR ]; then \
		mkdir $$DIR ; \
	fi ; \
	for file in `find 03_markdup-bowtie2-raw -type f -name '*.bam'`; do \
		fname=`basename $$file` ; \
		fname=$${fname%??????????} ; \
		for i in $$limits; do \
			qsub -b y -l h_vmem=8G -V -j y -m e -cwd -q all.q -N job-$$fname-rs$$i-markdup-bowtie2-raw "repeat-soaker -r $$GDIR -p $$i -o $$DIR/$$fname.$$i.bam $$file" ; \
		done ; \
    done

04_rs-markdup-bowtie2-tmatic:	03_markdup-bowtie2-tmatic
	LEN=19 ; \
	GDIR="/Volumes/hts_core/Shared/repeat_soaker/rmsk.hg19.bed" ; \
	DIR="04_rs-markdup-bowtie2-tmatic" ; \
	limits="0.75 0.50 0.25 0.01 0.00" ; \
	if [ ! -d $$DIR ]; then \
		mkdir $$DIR ; \
	fi ; \
	for file in `find 03_markdup-bowtie2-tmatic -type f -name '*.bam'`; do \
		fname=`basename $$file` ; \
		fname=$${fname%??????????} ; \
		for i in $$limits; do \
			qsub -b y -l h_vmem=8G -V -j y -m e -cwd -q all.q -N job-$$fname-rs$$i-markdup-bowtie2-tmatic "repeat-soaker -r $$GDIR -p $$i -o $$DIR/$$fname.$$i.bam $$file" ; \
		done ; \
    done
