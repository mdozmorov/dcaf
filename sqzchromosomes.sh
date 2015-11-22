n=$((`find . -type f -name "*.bed" | wc -l`)); \
i=0; \
for file in `find . -type f -name "*.bed"`; do \
	f=`basename $file`; \
	d=`dirname $file`; i=$((i+1)); \
	echo "Processing" $i "out of" $n ":" $file; \
	cat $file | grep "\bchr[0-9XYM][^_]\b" | \
		awk 'BEGIN {OFS="\t"} { if ( $3 <= $2) { print $1, $2, $2+1, $4, $5, $6 } else { print $1, $2, $3, $4, $5, $6 } }' | \
			sort -k1,1 -k2,2n -k3,3n | uniq > tmp.bed && mv tmp.bed $d/$f; \
	bgzip $file && tabix $file.gz; \
done

