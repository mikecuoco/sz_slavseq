import pysam
import subprocess

def prio_pair_rmdup(bam_in, txt):
    bam = pysam.AlignmentFile(bam_in, "rb") # read input bam file
    outfile = open(txt, "w") # open output file for writing

    r1 = None
    r2 = None

    while r2 = bam.fetch(until_eof = True):
        
        # filter out secondary matches (bwa mem compatibility)
	    # secondary matches won't be used in selection
	    # 256 = not primary alignment
	    # 2048 = supplementary alignment
	    # see https://broadinstitute.github.io/picard/explain-flags.html for SAM flags info
        if (r2.flag & 256) or (r2.flag & 2048):
            continue # continue to next iteration if true
        
        else if r1 != None:
            if (r2.qname != r1.qname) or ((r1.flag & 4) and (r1.flag & 64)):
                # if read names are not the same or 
                # if R1 is unmapped (4 = unmapped, 64 = first in pair), 
                # do nothing
            else: # this is a good pair
                if r1.flag & 128: # 128 = second in pair
                    # switch values of r1 and r2
                    tmp = r1
                    r1 = r2
                    r2 = tmp

                sumqual = 0
                # quality scores are at each base, so sum across bases of the read
                for i in range(0, r1.query_length):
                    sumqual += (r1.query_qualities[i] + r2.query_qualities[i])
                
                r1pos = r1.query_alignment_start + 1 # convert from 0-based to 1-based coordinates
                
                # label strand
                if r1.flag & 16:
                    r1strand = "-"
                else:
                    r1strand = "+"

                # join target ID, position, and strand info
                pos = ":".join([r1.reference_name, r1pos, r1strand])
                # join and print more metrics to output file
                all_fields = "\t".join([r1.qname, r1.mapping_quality + r2.mapping_quality, sumqual, pos]) + "\n"

                p1 = subprocess.Popen(["sort", "-S", "5000M", "-k", "4,4", "-k", "2,2rn", "-k", "3,3rn"], stdout=subprocess.PIPE)
                p1_out = p1.communicate(all_fields.encode())[0]
                p2 = subprocess.Popen(["uniq", "-f", "3", "-c"], stdin=p1_out, stdout=subprocess.PIPE)
                p1.stdout.close()

                p3 = subprocess.Popen(["perl", "-pe", 's/^ +(\\d+) +(\\S+)/\$2\\tXD:i:\$1/'], stdin=p2.stdout, stdout=subprocess.PIPE)
                p2.stdout.close()

                p4 = subprocess.Popen(["cut", "f", "1,2"], stdin=p3.stdout, stdout=subprocess.PIPE)
                p3.stdout.close()

                p5 = subprocess.run(["sort", "S", "5000M"], stdin=p4.stdout, stdout=outfile, text = True)
                p4.stdout.close()

                r2 = None
        
        r1 = r2

    bam.close()

prio_pair_rmdup("input.bam", "selected.txt")