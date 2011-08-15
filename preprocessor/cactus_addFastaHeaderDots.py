#!/usr/bin/env python

## Sequence names must be in "name.x" format in order for MAFJoin
## to work.  This script adds .0 to all sequence names in a FASTA
## file if they are not already in this format

import os
from optparse import OptionParser

from sonLib.bioio import fastaRead
from sonLib.bioio import fastaWrite
from sonLib.bioio import getTempFile

def fixHeader(header, num = "0"):
    suf = ""
    pref = header
    if '|1|' in header:
        idx = header.find('|1|')
        pref = header[:idx]
        suf = header[idx:]
    
    if '.' in pref:
        dotIdx = pref.find('.')
        if dotIdx == len(pref) - 1:
            pref = pref + num
    else:
        pref = pref + '.' + num
    
    return pref + suf
    
def main():
    ##########################################
    #Construct the arguments.
    ##########################################    
    
    usage = "usage: %prog [options] <fasta input file> <fasta output file>\n\n" + \
            "    <fasta file>:  fasta sequence to annotate\n"
    description = "Ensure sequence names are in name.x format by\n" + \
                    "adding .0's as necessary"
    parser = OptionParser(usage=usage, description=description)
    
    parser.add_option("--suffix", dest="suffix", type="string",
                      help="string to append after dot (default=0)",
                      default="0")
    
    options, args = parser.parse_args()
    
    if len(args) != 2:
        parser.print_help()
        return 1
    
    inputName = args[0]
    inputFile = open(inputName, "r")
    outputName = args[1]
    outputFile = open(outputName, "w")
     
    for header, seq in fastaRead(inputFile):
        fastaWrite(outputFile, fixHeader(header, options.suffix), seq)
            
    outputFile.close()
    inputFile.close()
    return 0
    
if __name__ == '__main__':
    exit(main())
