#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -M mdozmorov@vcu.edu
#PBS -N btbuild
#PBS -j oe
# PBS -o /home/glasser/Dissertation/Dozmorov/RNA-seq/cuffdiffout3

cd $PBS_O_WORKDIR

# Run this file from /home/mdozmorov/sequencing/data/ExtData/UCSC/hg38ext

FILE=chr1.fa,chr2.fa,chr3.fa,chr4.fa,chr5.fa,chr6.fa,chr7.fa,chr8.fa,chr9.fa,chr10.fa,chr11.fa,chr12.fa,chr13.fa,chr14.fa,chr15.fa,chr16.fa,chr17.fa,chr18.fa,chr19.fa,chr20.fa,chr21.fa,chr22.fa,chrM.fa,chrX.fa,chrY.fa,chrEBV.fa,chr1_KI270706v1_random.fa,chr1_KI270707v1_random.fa,chr1_KI270708v1_random.fa,chr1_KI270709v1_random.fa,chr1_KI270710v1_random.fa,chr1_KI270711v1_random.fa,chr1_KI270712v1_random.fa,chr1_KI270713v1_random.fa,chr1_KI270714v1_random.fa,chr2_KI270715v1_random.fa,chr2_KI270716v1_random.fa,chr3_GL000221v1_random.fa,chr4_GL000008v2_random.fa,chr5_GL000208v1_random.fa,chr9_KI270717v1_random.fa,chr9_KI270718v1_random.fa,chr9_KI270719v1_random.fa,chr9_KI270720v1_random.fa,chr11_KI270721v1_random.fa,chr14_GL000009v2_random.fa,chr14_GL000194v1_random.fa,chr14_GL000225v1_random.fa,chr14_KI270722v1_random.fa,chr14_KI270723v1_random.fa,chr14_KI270724v1_random.fa,chr14_KI270725v1_random.fa,chr14_KI270726v1_random.fa,chr15_KI270727v1_random.fa,chr16_KI270728v1_random.fa,chr17_GL000205v2_random.fa,chr17_KI270729v1_random.fa,chr17_KI270730v1_random.fa,chr22_KI270731v1_random.fa,chr22_KI270732v1_random.fa,chr22_KI270733v1_random.fa,chr22_KI270734v1_random.fa,chr22_KI270735v1_random.fa,chr22_KI270736v1_random.fa,chr22_KI270737v1_random.fa,chr22_KI270738v1_random.fa,chr22_KI270739v1_random.fa,chrUn_GL000195v1.fa,chrUn_GL000213v1.fa,chrUn_GL000214v1.fa,chrUn_GL000216v2.fa,chrUn_GL000218v1.fa,chrUn_GL000219v1.fa,chrUn_GL000220v1.fa,chrUn_GL000224v1.fa,chrUn_GL000226v1.fa,chrUn_KI270302v1.fa,chrUn_KI270303v1.fa,chrUn_KI270304v1.fa,chrUn_KI270305v1.fa,chrUn_KI270310v1.fa,chrUn_KI270311v1.fa,chrUn_KI270312v1.fa,chrUn_KI270315v1.fa,chrUn_KI270316v1.fa,chrUn_KI270317v1.fa,chrUn_KI270320v1.fa,chrUn_KI270322v1.fa,chrUn_KI270329v1.fa,chrUn_KI270330v1.fa,chrUn_KI270333v1.fa,chrUn_KI270334v1.fa,chrUn_KI270335v1.fa,chrUn_KI270336v1.fa,chrUn_KI270337v1.fa,chrUn_KI270338v1.fa,chrUn_KI270340v1.fa,chrUn_KI270362v1.fa,chrUn_KI270363v1.fa,chrUn_KI270364v1.fa,chrUn_KI270366v1.fa,chrUn_KI270371v1.fa,chrUn_KI270372v1.fa,chrUn_KI270373v1.fa,chrUn_KI270374v1.fa,chrUn_KI270375v1.fa,chrUn_KI270376v1.fa,chrUn_KI270378v1.fa,chrUn_KI270379v1.fa,chrUn_KI270381v1.fa,chrUn_KI270382v1.fa,chrUn_KI270383v1.fa,chrUn_KI270384v1.fa,chrUn_KI270385v1.fa,chrUn_KI270386v1.fa,chrUn_KI270387v1.fa,chrUn_KI270388v1.fa,chrUn_KI270389v1.fa,chrUn_KI270390v1.fa,chrUn_KI270391v1.fa,chrUn_KI270392v1.fa,chrUn_KI270393v1.fa,chrUn_KI270394v1.fa,chrUn_KI270395v1.fa,chrUn_KI270396v1.fa,chrUn_KI270411v1.fa,chrUn_KI270412v1.fa,chrUn_KI270414v1.fa,chrUn_KI270417v1.fa,chrUn_KI270418v1.fa,chrUn_KI270419v1.fa,chrUn_KI270420v1.fa,chrUn_KI270422v1.fa,chrUn_KI270423v1.fa,chrUn_KI270424v1.fa,chrUn_KI270425v1.fa,chrUn_KI270429v1.fa,chrUn_KI270435v1.fa,chrUn_KI270438v1.fa,chrUn_KI270442v1.fa,chrUn_KI270448v1.fa,chrUn_KI270465v1.fa,chrUn_KI270466v1.fa,chrUn_KI270467v1.fa,chrUn_KI270468v1.fa,chrUn_KI270507v1.fa,chrUn_KI270508v1.fa,chrUn_KI270509v1.fa,chrUn_KI270510v1.fa,chrUn_KI270511v1.fa,chrUn_KI270512v1.fa,chrUn_KI270515v1.fa,chrUn_KI270516v1.fa,chrUn_KI270517v1.fa,chrUn_KI270518v1.fa,chrUn_KI270519v1.fa,chrUn_KI270521v1.fa,chrUn_KI270522v1.fa,chrUn_KI270528v1.fa,chrUn_KI270529v1.fa,chrUn_KI270530v1.fa,chrUn_KI270538v1.fa,chrUn_KI270539v1.fa,chrUn_KI270544v1.fa,chrUn_KI270548v1.fa,chrUn_KI270579v1.fa,chrUn_KI270580v1.fa,chrUn_KI270581v1.fa,chrUn_KI270582v1.fa,chrUn_KI270583v1.fa,chrUn_KI270584v1.fa,chrUn_KI270587v1.fa,chrUn_KI270588v1.fa,chrUn_KI270589v1.fa,chrUn_KI270590v1.fa,chrUn_KI270591v1.fa,chrUn_KI270593v1.fa,chrUn_KI270741v1.fa,chrUn_KI270742v1.fa,chrUn_KI270743v1.fa,chrUn_KI270744v1.fa,chrUn_KI270745v1.fa,chrUn_KI270746v1.fa,chrUn_KI270747v1.fa,chrUn_KI270748v1.fa,chrUn_KI270749v1.fa,chrUn_KI270750v1.fa,chrUn_KI270751v1.fa,chrUn_KI270752v1.fa,chrUn_KI270753v1.fa,chrUn_KI270754v1.fa,chrUn_KI270755v1.fa,chrUn_KI270756v1.fa,chrUn_KI270757v1.fa,chrY_KI270740v1_random.fa

bowtie-build $FILE hg38ext
bowtie2-build $FILE hg38ext