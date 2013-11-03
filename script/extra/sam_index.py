#!/usr/bin/env python
import os, subprocess, glob
# files = os.listdir("*.bam") 
files = glob.glob('*.bam')
for file in files:
	# print file
	if os.path.exists(file+".bai"):
		continue
	try:
		cmd = "samtools index "+file
		subprocess.call(cmd, shell=True)
		print "Processed "+file
	except Exception,e:
		print "Error "+e
