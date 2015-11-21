#!/bin/bash

cd "$( dirname "${BASH_SOURCE[0]}" )"

# be sure to start with a fresh & known enviroment (hopefully)
. /usr/local/Modules/default/init/bash
module unload bwa
module unload samtools
module unload picard
module unload gatk
module unload fastqc
module unload tabix
module unload btags
module load bwa/0.7.10
module load samtools/0.1.19
module load picard/1.99
module load gatk/3.2-2
module load fastqc/0.10.1
module load tabix/0.2.6
module load btangs/1.6.0

set -o pipefail


# It is easier to use variables for thse paths
GATK_REF=/Volumes/hts_core/Shared/gatk_resources/2.3/b37/bwa_7_6_indexed/human_g1k_v37.fasta
GATK_DBSNP=/Volumes/hts_core/Shared/gatk_resources/2.3/b37/dbsnp_137.b37.vcf
GATK_BIN=`which gatk`
GATK_BASE=`dirname ${GATK_BIN}`"/.."

SAMPLE="test"

FASTQ1=/net/flotsam.nfs.ngs.omrf.in/ifs/groups/wren_lab/scratch/temp/ERR034546_1.filt.fastq
FASTQ2=/net/flotsam.nfs.ngs.omrf.in/ifs/groups/wren_lab/scratch/temp/ERR034546_2.filt.fastq

function final_clean() {
rm -rf "$TMP_DIR"
}
TMP_DIR=/tmp
if [ "$?" -ne "0" ]; then
  echo "Failed to make our temp work dir"
  exit 1
fi
trap final_clean EXIT
export GATK_JAVA_OPTS="-Djava.io.tmpdir=${TMP_DIR}"

if [ ! -e 05_dup_marked/cleaned.bam ]; then
# initial PCR clean by 'bins' via btangs


# fastqc info
#mkdir qc
#qsub  -p -1000 -l virtual_free=2G,h_vmem=4G -o logs -b y -V -j y -cwd -N a_test_qc fastqc -o qc ${FASTQ1} ${FASTQ2}

# setup input sams, will get illumina scores to standard sanger
export JAVA_MEM_OPTS="-Xmx16G"

# Now take those & actually generate the alignment SAM output, paired or single
mkdir 03_sorted_bams
qsub  -pe threaded 12 -l virtual_free=1G,mem_free=1G,h_vmem=48G -o logs -sync y -t 1-1 -b y -V -j y -cwd -N a_test_bwa_alignment bwa_mem_qsub_tasked.rb ${TMP_DIR} 03_sorted_bams /Volumes/hts_core/Shared/gatk_resources/2.3/b37/bwa_7_6_indexed/human_g1k_v37.fasta paired '"@RG\\tID:test_001_s_1\\tSM:test\\tPL:Illumina\\tPU:1"' ${FASTQ1} ${FASTQ2}

if [ "$?" -ne "0" ]; then
  echo -e "Failure with bwa alignment"
  exit 1
fi


# Now we might have had many input SAMs, so let us merge those all into a single BAM using picard
# While we do that we shall also sort it, mark possible duplicates & make an index
mkdir 05_dup_marked
qsub  -l virtual_free=8G,mem_free=8G,h_vmem=56G -o logs -sync y -b y -V -j y -cwd -N a_test_merge_mark_dups \
  picard MarkDuplicates TMP_DIR=${TMP_DIR} INPUT=03_sorted_bams/0.bam  \
  OUTPUT=./05_dup_marked/cleaned.bam METRICS_FILE=./05_dup_marked/mark_dups_metrics.txt \
  VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=6000000 CREATE_INDEX=True COMPRESSION_LEVEL=8

if [ "$?" -ne "0" ]; then
  echo -e "Failure with marking the duplicates"
  exit 1
fi


rm -rf 03_sorted_bams


# TODO stop here if pre_gatk only
if [ "$PRE_GATK_ONLY" == "Y" ]; then
rm -rf 00_inputs \
01_bwa_aln_sai \
02_bwa_alignment 03_sorted_bams \
${TMP_DIR}
echo noop

rm -f qc/*.zip

touch finished.txt
exit 0
fi
else
  rm -f finished.txt
  TMP_DIR=/tmp
if [ "$?" -ne "0" ]; then
  echo "Failed to make our temp work dir"
  exit 1
fi
fi #if 05_dup_marked/cleaned.bam already existed
# start here if pre_gatk already done

    # Calculate intervals for realignment
    mkdir 06_intervals
    qsub  -pe threaded 6 -R y -o logs -sync y -b y -V -j y -cwd -N a_test_intervals \
     -l virtual_free=1G,mem_free=1G,h_vmem=20G \
     gatk -T RealignerTargetCreator -R ${GATK_REF} -I ./05_dup_marked/cleaned.bam -o ./06_intervals/cleaned.intervals -nt 10

    if [ "$?" -ne "0" ]; then
     echo -e "Failure with target realigment creation"
     exit 1
    fi

    # Now realign & fix any mate info
    mkdir 07_realigned_bam
    unset JAVA_MEM_OPTS
    qsub  -o logs -sync y -b y -V -j y -cwd -N a_test_realign \
     -l virtual_free=5G,mem_free=4G,h_vmem=8G \
     gatk -T IndelRealigner -known /Volumes/hts_core/Shared/gatk_resources/2.3/b37/1000G_phase1.indels.b37.vcf -known /Volumes/hts_core/Shared/gatk_resources/2.3/b37/Mills_and_1000G_gold_standard.indels.b37.vcf -R ${GATK_REF} -I ./05_dup_marked/cleaned.bam \
     --targetIntervals ./06_intervals/cleaned.intervals -o ./07_realigned_bam/cleaned.bam --maxReadsInMemory 1000000 

    if [ "$?" -ne "0" ]; then
     echo -e "Failure with indel realigmnent"
     exit 1
    fi


  # BaseRecalibrator
  mkdir 08_uncalibated_covariates
  unset JAVA_MEM_OPTS
  qsub  -pe threaded 6 -R y -o logs -sync y -b y -V -j y -cwd -N a_test_bqsr \
   -l virtual_free=1G,mem_free=4G,h_vmem=8G \
   gatk -T BaseRecalibrator -R ${GATK_REF} --knownSites /Volumes/hts_core/Shared/gatk_resources/2.3/b37/dbsnp_137.b37.vcf --knownSites /Volumes/hts_core/Shared/gatk_resources/2.3/b37/Mills_and_1000G_gold_standard.indels.b37.vcf -I ./07_realigned_bam/cleaned.bam \
   -o ./08_uncalibated_covariates/recal_data.grp -nct 6

  if [ "$?" -ne "0" ]; then
   echo -e "Failure counting covariates"
   exit 1
  fi

  mkdir 10_recalibrated_bam
  unset JAVA_MEM_OPTS
  qsub  -pe threaded 6 -R y -o logs -sync y -b y -V -j y -cwd -N a_test_recalibrate \
   -l virtual_free=1G,mem_free=4G,h_vmem=8G \
   gatk -T PrintReads -R ${GATK_REF} -I ./07_realigned_bam/cleaned.bam -BQSR ./08_uncalibated_covariates/recal_data.grp \
   -o ./10_recalibrated_bam/recalibrated.bam --bam_compression 8 -nct 6

  if [ "$?" -ne "0" ]; then
   echo -e "Failure reclibrating bam"
   exit 1
  fi

  mkdir 13_final_bam

  mv ./10_recalibrated_bam/recalibrated.bam ./13_final_bam/test.bam
  mv ./10_recalibrated_bam/recalibrated.bai ./13_final_bam/test.bam.bai


# fastqc info
#qsub  -p -1000 -l virtual_free=2G,h_vmem=4G -o logs -b y -V -j y -cwd -N a_test_qc fastqc -o qc ./13_final_bam/test.bam


# Clean up after ourselves
rm -rf 00_inputs \
01_bwa_aln_sai \
02_bwa_alignment 03_sorted_bams \
05_dup_marked \
06_intervals \
07_realigned_bam \
08_uncalibated_covariates \
10_recalibrated_bam \
11_calibated_covariates



# Get some summary stats on the alignment off in the background
JAVA_MEM_OPTS="-Xmx4G" qsub  -l virtual_free=2G,mem_free=2G,h_vmem=6G -o logs -b y -V -j y -cwd -N a_test_alignment_summary \
picard CollectAlignmentSummaryMetrics INPUT=13_final_bam/test.bam OUTPUT=13_final_bam/align_summary.txt VALIDATION_STRINGENCY=LENIENT

    # Finally Haplotypecaller in gVCF mode or is Gvcf mode

    export JAVA_MEM_OPTS="-Xmx24G"
    qsub  -pe threaded 4 -o logs -sync y -b y -V -j y -cwd -N a_test_variants \
    -l virtual_free=3G,mem_free=3G,h_vmem=28G gatk -T HaplotypeCaller \
    --pair_hmm_implementation VECTOR_LOGLESS_CACHING -ERC GVCF -nct 4 -R ${GATK_REF} \
    -I ./13_final_bam/test.bam -o test.gvcf \
    -variant_index_type LINEAR -variant_index_parameter 128000 -D ${GATK_DBSNP} -L /net/flotsam.nfs.ngs.omrf.in/ifs/groups/wren_lab/scratch/temp/refGene.bed
    # -stand_emit_conf 10.0 -stand_call_conf 30.0

    if [ "$?" -ne "0" ]; then
     echo -e "Failure GVCF"
     exit 1
    fi

  bgzip test.gvcf
  tabix -p vcf test.gvcf.gz




# gzip & clean up some of the cleaned input as we keep some it for now, but don't need it fullsized
if [ "$PRE_GATK_ONLY" != "Y" ]; then
echo noop
fi

rm -f qc/*.zip
rm -rf ${TMP_DIR}

touch finished.txt
