# A collection of miscellaneous scripts
Will gradually be documented or tossed

alignDirectory
annotationAnalysis
appendUniquePostfix.py

`bedgroomer.py` - Check if the end coordinate in a bed file is at least 1bp larger, and fixes if not. Useful to groom files with SNP's coordinates, as they may contain deletions with equal start and end coordinates. Can be done with awk.

bigWigCount
bigWigDescribe
categorize.py
collapseRows

`convert.py` - Reads 3-column file containing "line #", "Gene(s)", "Ratio"; Genes may be " /// " separated; Splits gene line into a set of genes; Outputs "Gene" - "Ratio" pair. For Affy data

count-bam
exome2.py
exportSOFT
extract.sh

`extract_UCSC.py` - Connects to MySQL database, extracts BroadHMM subgroups into separate files based on 'name' column. Does not work on python3 - MySQL module not supporrted.

fetchGOMapping
fetchGeneLoci

`find_dupes.sh` - Finds duplicate files and creates a script to remove them

ftp.sh
ftp_download.py
generate_random_bed.sh
gff2bed.py
group.py
grtk.sh

`gwas2bed.py` - Extract genomic coordinates of individual traits from GWAScatalog, and output them into trait-specific .bed files.

gwas_subsets.py
gwassort.sh
hgExport
hypergeom.py
immgen.py
import-soft
lineCount
mat2vw
merge.py
mergeT.py
mysql2sqlite.sh
ncbi-gff-to-bed
ngram.py
nonzero.py
orphanfinder.py
overlapStatistics
pattern_scan.py
prepareModel
promoter_extract
rand.py
rand.txt
rdover.py

`repeatSoakerParser.py` - Read in FASTA file and concatenate FASTA sequences interspersed by randomly generated sequences of fixed length.

`remspaces.sh` - remove spaces from file names (renaming)

rfam2bed.py
rm_missing.py
rmsemi.py
sam_index.py

`sqzchromosomes.sh` - filters out non-canonical chromosomes from BED files

`sqzsinglequotes.sh` - remove single quotes from file names (renaming)

`sqzspaces.sh` - remove spaces from file names (renaming)

toBED
toBigWig

`tomatrix.py` - convert .ped and .map PLINK files to a matrix patients x SNPs

transpose
transpose.py
transpose2.py
vcf_to_gene.sh
vcfparser.py
vwCrossValidate
vwFillTemplate
xyz.py