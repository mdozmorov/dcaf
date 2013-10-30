#!/usr/bin/env python
import sys, os, subprocess
with open(sys.argv[1]) as accessions_file:
	accessions = set(line.strip() for line in accessions_file.readlines())
target_folder = sys.argv[2]
if not os.path.exists(target_folder):
	os.makedirs(target_folder)

from ftplib import FTP
ftp = FTP('ftp-trace.ncbi.nih.gov')
ftp.login()
ftp.cwd("1000genomes/ftp/data")

BASE_DIR="ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data"
CHROM_LOC="12:10997448-11549498"
for accession in accessions:
	print accession
	# target_file = "%s/%s.PRH1_exome.bam" % (target_folder, accession)
	target_file = "%s/%s.PRH1_alignment.bam" % (target_folder, accession)
	if os.path.exists(target_file):
		continue
	try:
		# files = ftp.nlst(accession+"/exome_alignment")
		files = ftp.nlst(accession+"/alignment")
		file = filter(lambda f: ".mapped." in f and f.endswith(".bam"), files)[0]
		# file = [f for f in files if ".mapped." in f and f.endswith(".bam")][0]
		cmd = "samtools view -bh %s/%s %s > %s" % (BASE_DIR, file, CHROM_LOC, target_file)
		subprocess.call(cmd, shell=True)
	
	except Exception:
		print "Oops! Something wrong with sample "+accession
ftp.quit()
