#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

#nknguyen@soe.ucsc.edu
#09/24/2010
#
#Get the reference sequence from maf output of cactus_MAFGenerator (with 
#options --orderByReference and --includeReferenceSequence specified)
#Input: cactus.maf Output refseq.fa
#Given that blocks are in order of the 'reference' in cactus.maf, and 'reference'
#is always positive.

import os
import sys
import re
from optparse import OptionParser

def getRefSeq():
    seq = ''
    #count = 0
    #totalbases = 0
    #go through input file (stdin) and get the sequences
    for line in sys.stdin.readlines():
        p = re.compile('^s\treference')
        if p.match(line):
            #0   1      2      3       4         5          6
            #s, name, start, length, strand, totalLength, sequence
            list = line.rstrip().split('\t')
	    if len(list) != 7:
	        sys.stderr.write("Line %s has wrong format\n" %line)
	        continue
            else:
                seq += list[6]
                #count += 1
                #totalbases += int(list[3])
    sys.stdout.write(">reference\n")
    #sys.stderr.write("Count: %d; totalBase: %d; Total bases in the ref: %d\n" %(count, totalbases, len(seq)))
    for i in range(0, len(seq), 50):
        sys.stdout.write(seq[i:i+50] + '\n')

if __name__ == "__main__" :
    usage = "usage: %prog [options] < input maf file > output fasta file"
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()
    getRefSeq()

