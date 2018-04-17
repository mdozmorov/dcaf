# Scripts for basic RNA-seq data processing and summarization

`subread-align` (http://bioinf.wehi.edu.au/subread/) and `featureCounts` pipeline (http://bioinf.wehi.edu.au/featureCounts/)

## Recommended folder structure

- `00_raw` - place gzipped FASTQ files here.
- `00_fastqc-raw` - FASTQC quality control output. Created by `submit00_fastqc.sh`.
- `01_trimmed` - adapter-trimmed FASTQ data.
- `02_subread-align` - BAM files aligned to the reference genome.
- `03_sorted` - sorted BAM files.
- `03_fastqc-bam` - FASTQC of aligned BAM files.

## Files

- `fastq-summary` - Python script to summarize multiple outputs from FASTQC into one tab-separated text file
- `get_sra_data.sh` - example of SRA data download
- `omicsoft.fa` - FASTA sequences of adapters
- `trimmomatic-0.33.jar` - make Trimmomatic jar file available from the local folder

## Scripts

- `submit00_fastqc.sh` - FASTQC quality control.
- `submit00_RseQC_geneBody_coverage.sh` - sequencing coverage across gene bodies, using `RseQC`
- `submit01_bowtie-build_hg38ext.sh`, `submit01_bowtie-build.sh` - examples of building Bowtie index
- `submit01_cutadapt.sh` - trim adapters using `Cutadapt`. Adjust for single- or paired-end reads.
- `submit01_trimmomatic.sh` - trim adapters using `Trimmomatic`. Adjust for single- or paired-end reads.
- `submit02_subindex.sh` - creates reference genome index.
- `submit02_subread.sh` - aligns trimmed FASTQ files to a reference genome. Adjust for single- or paired-end reads.
- `submit03_samsort.sh` - sorts aligned BAM files.
- `submit02_tophat.sh` - TopHat run
- `submit03_fastqc_bam.sh` - FASTQC quality control of aligned BAM files.
- `submit03_flagstat.sh` - flafstat stats of BAM files
- `submit03_picard.sh` - Picard tools run
- `submit03_featureCounts.sh` - summarize gene counts using `featureCounts`
- `submit03_unmapped.sh` - extract unmapped reads
- `submit04_cuffnorm.sh` - create expression matrix normalized to library size, last step after cufflinks
- `submit04_DEXSeq_prepare_annotation.sh` - prepare annotations for DEXseq analysis
- `submit04_DEXSeq_count.sh` - get exon counts using DEXSeq annotation
- `submit04_htseq.sh` - summarize gene counts using `htseq-count`
