Miscellaneous scripts
---------------------
Will gradually be documented or tossed

`btang.sh` - duplicate FASTQ remover

`fastqc-summary` - summarizes multipe *.zip FASTQC outputs into text matrix

`join_pairs.sh` - duplicate FASTQ remover

`makemtx.py` - converts cuffdiff output into gene expression matrix

omicsoft.fa
README.md
submit_bam2wig.sh
submit_cutadapt.sh
submit_star.sh
submit00_btangs.sh
submit00_fastqc.sh
submit00_raw_merge.sh
submit01_fastqc.sh
submit01_trimmomatic.sh
submit02_fastqc.sh
submit02_sickle.sh
submit03_bowtie-sickle.sh
submit03_bowtie-tmatic.sh
submit03_picard.sh
submit03_rna-seqc.sh
submit03_samtools.sh
submit03_star_rsoak.sh
submit04_cufflinks.sh
submit04_htseq.sh
submit05_homer.sh
submit05_macs2-bowtie.sh
submit05_samtools.sh
submit05_sicer.sh
submit06_multicov_macs2.sh
submit08_samtools_dup.sh
submit1.sh
submit1.stranded.sh
submit1.unstranded.sh
trimends

## Getting uniquely mappable regions from Umap

```
# Download manually: https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv
wget https://www.pmgenomics.ca/hoffmanlab/proj/bismap/raw/hg38/k36.umap.multitrackmappability.wg.gz
wigToBigWig -clip k36.umap.multitrackmappability.wg GRCh38_EBV.chrom.sizes.tsv k36.umap.multitrackmappability.bw
bigWigToBedGraph k36.umap.multitrackmappability.bw k36.umap.multitrackmappability.bedGraph
awk 'BEGIN {FS=OFS="\t"} {if ($4 >= .75) print $1, $2, $3}' k36.Umap.MultiTrackMappability.bedGraph | bedtools merge > k36.Umap.MultiTrackMappability.filtered.bed
```

Produces a BED file qith ~19M entries, ~80bp long.