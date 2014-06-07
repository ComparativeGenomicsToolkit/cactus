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
    
    parser.add_option("--checkAlphaNumeric", dest="checkAlphaNumeric", action="store_true",
                      help="Checks that the first word contains only alphanumeric characters, periods or underscores.",
                      default=False)

    parser.add_option("--checkUCSC", dest="checkUCSC", action="store_true",
                      help="Checks that suffix of the first word after the last '.' character contains only alpha-numeric characters or underscores and is unique.",
                      default=False)

    options, args = parser.parse_args()
    
    if len(args) != 1:
        parser.print_help()
        return 1
    
    inputName = args[0]
    inputFile = open(inputName, "r")
     
    seen = set()
    for header, seq in fastaRead(inputFile):
        mungedHeader = header.split()[0]
        if options.checkAlphaNumeric and "".join([ i for i in mungedHeader if str.isalnum(i) ]) != mungedHeader: #Check is only alpha numeric
            raise RuntimeError("We found a non-alpha numeric character in the fasta header, and the config file (checkAlphaNumeric option) demands that all fasta headers be alpha numeric: %s" % header)
        if options.checkUCSC:
            mungedHeader = mungedHeader.split('.')[-1]
            if "".join([ i for i in mungedHeader if (str.isalnum(i) or i == '_' or i == '-' or i == ':') ]) != mungedHeader:
                raise RuntimeError("We found a non-alpha numeric, '-', ':' or '_' prefix in the fasta header (UCSC Names option), if you want to be able to generate a UCSC assembly hub to display your alignment then please modify the first word after the '>' and before the first '.' in every fasta header to only contain alpha-numeric, '_', ':' or '-' characters. The offending header: %s" % header)
        if mungedHeader in seen:
            raise RuntimeError("We found a duplicated fasta header, the first word of each fasta header should be unique within each genome, as this is a requirement for the output HAL file or any MAF file subsequently created. Please modify the input fasta file. Offending duplicate header: %s" % header)
        seen.add(mungedHeader)
    inputFile.close()
    return 0
    
if __name__ == '__main__':
    exit(main())