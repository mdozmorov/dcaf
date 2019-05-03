# Input folder
DIRIN=03_sorted-star/

# Output folder
DIROUT=04_MACS2-star

mkdir -p $DIROUT


# Output file
OUTPUT="PF20.bed"

macs2 callpeak -t $DIRIN"anti_PF20_0319_S57_L006_R1_001.Aligned.out.bam" $DIRIN"anti_PF20_0320_S60_L006_R1_001.Aligned.out.bam" $DIRIN"anti_PF20_0326_S63_L006_R1_001.Aligned.out.bam" -c $DIRIN"anti_IgG_0319_S55_L006_R1_001.Aligned.out.bam" $DIRIN"anti_IgG_0320_S58_L006_R1_001.Aligned.out.bam" $DIRIN"anti_IgG_0326_S61_L006_R1_001.Aligned.out.bam" -g mm --outdir $DIROUT --name $OUTPUT --verbose 3

# Output file
OUTPUT="SOX5.bed"

macs2 callpeak -t $DIRIN"anti_SOX5_0319_S56_L006_R1_001.Aligned.out.bam" $DIRIN"anti_SOX5_0320_S59_L006_R1_001.Aligned.out.bam" $DIRIN"anti_SOX5_0326_S62_L006_R1_001.Aligned.out.bam" -c $DIRIN"anti_IgG_0319_S55_L006_R1_001.Aligned.out.bam" $DIRIN"anti_IgG_0320_S58_L006_R1_001.Aligned.out.bam" $DIRIN"anti_IgG_0326_S61_L006_R1_001.Aligned.out.bam" -g mm --outdir $DIROUT --name $OUTPUT --verbose 3

