#!/usr/bin/env python

## Filter small sequences out of a fasta file.  For use with flies,
## for example, where scaffolds of length <200kb seem to be considered
## no mans land

import os
from optparse import OptionParser

from sonLib.bioio import fastaRead
from sonLib.bioio import fastaWrite
from sonLib.bioio import getTempFile

def doFilter(header, seq, options):
    return header.find(options.prefix) == 0 and len(seq) < options.length
    
def main():
    ##########################################
    #Construct the arguments.
    ##########################################    
    
    usage = "usage: %prog [options] <fasta input file> <fasta output file>\n\n" + \
            "    <fasta file>:  fasta sequence to filter\n"
    description = "Ensure sequences have length >= length\n"
                    
    parser = OptionParser(usage=usage, description=description)
    
    parser.add_option("--prefix", dest="prefix", type="string",
                      help="only filter sequences with prefix in name",
                      default="")
    parser.add_option("--length", dest="length", type="int",
                      help="filter shorter than length [default=1000]",
                      default=1000)
    
    options, args = parser.parse_args()
    
    if len(args) != 2:
        parser.print_help()
        return 1
    
    inputName = args[0]
    inputFile = open(inputName, "r")
    outputName = args[1]
    outputFile = open(outputName, "w")
     
  
    for header, seq in fastaRead(inputFile):
        if doFilter(header, seq, options) == False:
            fastaWrite(outputFile, header, seq)
      
    outputFile.close()
    inputFile.close()  
    return 0

if __name__ == '__main__':
    exit(main())
