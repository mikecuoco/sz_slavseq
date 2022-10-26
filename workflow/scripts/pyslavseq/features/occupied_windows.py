#!/usr/bin/env python
__author__ = "Apuã Paquola"

import subprocess
from ..genome import Interval


def occupied_windows_in_genome(genome, window_size, step, tabix_filename):
    occ = dict()
    cmd = "gunzip -c " + tabix_filename + "| cut -f 1-3,5,8 | awk 'and($4,64) && $5>=40 {print}' | cut -f 1-3 " #| head -5
    f = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    for line in f.stdout:
        a = line.decode().rstrip('\r\n').split(sep='\t')
        
        if a[0] not in occ:
            occ[a[0]]=set()
            
        wstart = max(-1, int(a[1]) - window_size)
        wend = int(a[2])
        
        for x in range(wstart // step + 1, (wend-1) // step + 1):
            occ[a[0]].add(x*step)
        
    for chrom in sorted(occ.keys()):
        for pos in sorted(occ[chrom]):
            yield genome.fit_interval(Interval(chrom,
                                               pos,
                                               pos + window_size))
