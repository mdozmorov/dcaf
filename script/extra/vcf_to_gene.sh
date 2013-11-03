#/bin/bash
#i=1
for file in output/*.bed
do
#	i=$(($i+1))
#	echo "111"> $file$i".txt"
	intersectBed -wa -wb -a $file -b exons2genes.bed | cut -f7 | sort | uniq -c | awk '{print $2,$1}' | sort -n -r -k2 > $file".rnk"
done
