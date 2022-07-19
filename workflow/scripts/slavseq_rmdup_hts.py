#!/usr/bin/env python

import pysam
import subprocess
import sys
import os
import locale

# def use_tmp_dir(path):
#     # TODO
#     tmpbase = path
#     old = path

def prio_pair_rmdup(filename, out_filename):
    bam = pysam.AlignmentFile(filename, "rb") # read input bam file
    outfile = open(out_filename, "a+") # open output file for writing

    r1 = None
    r2 = None

    for r2 in bam.fetch(until_eof = True):
        
        # filter out secondary (and supplementary) matches (bwa mem compatibility)
	    # secondary matches won't be used in selection
        if (r2.is_secondary) or (r2.is_supplementary):
            continue # continue to next iteration if true

        if r1 != None:
            if not (r2.qname != r1.qname or (r1.is_unmapped and r1.is_read1)):
                # if read names are not the same or if R1 is unmapped, do nothing
                # otherwise, we have a good pair

                if r1.is_read2:
                    # switch values of r1 and r2
                    tmp = r1
                    r1 = r2
                    r2 = tmp

                sumqual = 0
                # quality scores are at each base, so sum across bases of each read
                query_qualities = r1.query_qualities + r2.query_qualities
                for q in query_qualities:
                    sumqual += q
                
                r1pos = r1.query_alignment_start + 1 # convert from 0-based to 1-based coordinates
                
                # label strand
                r1strand = "-" if r1.is_reverse else "+"

                # join target ID, position, and strand info
                pos = ":".join([r1.reference_name, str(r1pos), r1strand])
                # join and print more metrics to output file
                all_fields = "\t".join([r1.qname, str(r1.mapping_quality + r2.mapping_quality), str(sumqual), pos])
                
                subprocess.run(
                    "sort -S 5000M -k 4,4 -k 2,2rn -k 3,3rn | uniq -f 3 -c |" +
                    "perl -pe 's/^ +(\\d+) +(\\S+)/\$2\\tXD:i:\$1/' | cut -f 1,2 | sort -S 5000M", 
                    input=all_fields, stdout=outfile, shell=True, check=True, text=True)

                r2 = None
        
        r1 = r2

    outfile.close()

input_bam_fn = sys.argv[1]
output_bam_fn = sys.argv[2]

# input_bam_path = os.path.abspath(sys.argv[1])

if os.path.exists(output_bam_fn):
    sys.exit("Output file already exists!")

locale.setlocale(locale.LC_ALL, "C")

# pwd = os.getcwd().strip()
# curdir, tmpdir = use_tmp_dir(pwd)
# os.symlink(input_bam_path, tmpdir + "/input.bam")
# os.chdir(tmpdir)

prio_pair_rmdup(
    input_bam_fn, 
    "results/selected.txt")

input_bam = open(input_bam_fn, "r")
header = open("results/header.txt", "w+")
output_bam = open(output_bam_fn, "w+")

subprocess.run(["samtools", "view", "-H"], stdin=input_bam, stdout=header, check=True, text=True)

input_bam.seek(0) # reset file pointer to start of file
p1 = subprocess.Popen(["samtools", "view"], stdin=input_bam, stdout=subprocess.PIPE, text=True)
p2 = subprocess.run(
    "sort -T ./ -S 1500M -s -k 1,1 |" +
    "join -t '\t' - results/selected.txt | cat results/header.txt - | samtools view -S -b -",
    stdin=p1.stdout, stdout=output_bam, shell=True, check=True, text=True)

input_bam.close()
header.close()
output_bam.close()

os.remove("results/selected.txt")
os.remove("results/header.txt")