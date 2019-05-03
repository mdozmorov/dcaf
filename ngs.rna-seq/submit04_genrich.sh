# DIRIN=03_sorted
# DIROUT=03_sorted.n
# mkdir -p $DIROUT
# # Sort by name
# for file in $DIRIN"/"*.bam; do sambamba sort -n -p $file -o $DIROUT"/"`basename $file`; done

# Input folder
DIRIN=03_sorted.n/

# Output folder
DIROUT=04_Genrich

mkdir -p $DIROUT


# Control file
CONTROL=$DIRIN"anti_IgG_0319_S55_L006_R1_001.bam",$DIRIN"anti_IgG_0320_S58_L006_R1_001.bam",$DIRIN"anti_IgG_0326_S61_L006_R1_001.bam"

# Test file
TEST=$DIRIN"anti_PF20_0319_S57_L006_R1_001.bam",$DIRIN"anti_PF20_0320_S60_L006_R1_001.bam",$DIRIN"anti_PF20_0326_S63_L006_R1_001.bam"
# Output file
OUTPUT=$DIROUT"/PF20.bed"

./Genrich -t $TEST -o $OUTPUT -c $CONTROL


# Test file
TEST=$DIRIN"anti_SOX5_0319_S56_L006_R1_001.bam",$DIRIN"anti_SOX5_0320_S59_L006_R1_001.bam",$DIRIN"anti_SOX5_0326_S62_L006_R1_001.bam"
# Output file
OUTPUT=$DIROUT"/SOX5.bed"

./Genrich -t $TEST -o $OUTPUT -c $CONTROL

