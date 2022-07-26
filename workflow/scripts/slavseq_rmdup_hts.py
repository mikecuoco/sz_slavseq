#!/usr/bin/env python
__author__ = "Rohini Gadde"

import pysam
import subprocess
import sys
import os
import argparse
import tempfile
import shutil
import glob
from pathlib import Path

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
                
                r1pos = r1.reference_end + 1 if r1.is_reverse else r1.reference_start + 1 # convert from 0-based to 1-based coordinates
                
                # label strand
                r1strand = "-" if r1.is_reverse else "+"

                # join target ID, position, and strand info
                pos = ":".join([r1.reference_name, str(r1pos), r1strand])
                
                # join and print more metrics to output file
                all_fields = "\t".join([r1.qname, str(r1.mapping_quality + r2.mapping_quality), str(sumqual), pos])
                
                outfile.write(all_fields + "\n")

                r2 = None
        
        r1 = r2

    outfile.close()
    
def parse_args():

    parser = argparse.ArgumentParser(description="Remove duplicates from BAM files")
    parser.add_argument(
        "-b", "--bam", type=Path, required=True, help="input BAM file"
    )
    parser.add_argument(
        "-o", "--out", type=Path, required=True, help="output BAM file"
    )
    
    args = parser.parse_args()

    return args

def main():
    
    # get arguments
    args = parse_args()

    input_bam_fn = args.bam
    output_bam_fn = args.out

    if os.path.exists(output_bam_fn):
        sys.exit("Output file already exists!")

    os.environ['LC_ALL'] = "C"

    curdir = os.getcwd()
    tmpdir = tempfile.mkdtemp(dir="./")
    os.symlink(input_bam_fn, tmpdir + "/input.bam")
    os.chdir(tmpdir)

    prio_pair_rmdup(
        "input.bam", 
        "all_fields.txt")
    
    unsorted_ids = open("all_fields.txt", "r")
    sorted_ids = open("selected.txt", "w+")

    # NOTE: shell=True can lead to security vulnerabilities
    subprocess.run(
        'sort -S 5000M -k 4,4 -k 2,2rn -k 3,3rn | uniq -f 3 -c | ' +
        'perl -pe "s/^ +(\\d+) +(\\S+)/\$2\\tXD:i:\$1/" | cut -f 1,2 | sort -S 5000M',
        stdin=unsorted_ids, stdout=sorted_ids, shell=True, check=True, text=True)
    
    input_bam = open("input.bam", "r")
    header = open("header.txt", "w+")
    output_bam = open("output.bam", "w+")

    subprocess.run(["samtools", "view", "-H"], stdin=input_bam, stdout=header, check=True)

    input_bam.seek(0) # reset file pointer to start of file
    p1 = subprocess.Popen(["samtools", "view"], stdin=input_bam, stdout=subprocess.PIPE)
    p2 = subprocess.run(
        'sort -T ./ -S 1500M -s -k 1,1 | ' +
        'join -t "\t" - selected.txt | cat header.txt - | samtools view -S -b -',
        stdin=p1.stdout, stdout=output_bam, shell=True, check=True)

    input_bam.close()
    header.close()
    output_bam.close()

    os.chdir(curdir)
    shutil.move(tmpdir + "/output.bam", output_bam_fn)

if __name__ == "__main__":
    main()
    
    tmpdirs = glob.glob("tmp*")
    if len(tmpdirs) != 0:
        for tmp in tmpdirs:
            shutil.rmtree(tmp)