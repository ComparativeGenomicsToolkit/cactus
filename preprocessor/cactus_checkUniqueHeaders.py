#!/usr/bin/env python

"""Checks headers are all unique.
"""

import os
from optparse import OptionParser

from sonLib.bioio import fastaRead
from sonLib.bioio import fastaWrite
from sonLib.bioio import getTempFile

def main():
    ##########################################
    #Construct the arguments.
    ##########################################    
    
    usage = "usage: %prog [options] <fasta input file>\n\n" + \
            "    <fasta file>:  fasta sequence to check for unique headers\n"
    description = "Ensure sequence names are unique\n" 
    parser = OptionParser(usage=usage, description=description)

    options, args = parser.parse_args()
    
    if len(args) != 1:
        parser.print_help()
        return 1
    
    inputName = args[0]
    inputFile = open(inputName, "r")
     
    seen = set()
    for header, seq in fastaRead(inputFile):
        if header in seen:
            raise RuntimeError("We found a duplicated fasta header: %s" % header)
        seen.add(header)
    inputFile.close()
    return 0
    
if __name__ == '__main__':
    exit(main())