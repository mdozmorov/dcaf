#!/bin/bash

for file in `ls gwascatalog_bed/*`; do
	echo $file
	# wc -l $file | cut -d " " -f1 -
	if [ ! `wc -l "$file" |cut -d " " -f1 -` -ge 50 ]
	then
		mv "$file" gwascatalog_bed_less50
	else
		mv "$file" gwascatalog_bed_more50
	fi
done

