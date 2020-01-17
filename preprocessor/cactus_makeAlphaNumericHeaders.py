#!/usr/bin/env python3

"""Removes all non-alpha-numeric chatavyers from the fasta headers of a fasta file.
"""

import os
from optparse import OptionParser

from sonLib.bioio import fastaRead
from sonLib.bioio import fastaWrite
from sonLib.bioio import getTempFile

def fixHeader(header):
    return "".join([ i for i in header if str.isalnum(i) ])

def main():
    ##########################################
    #Construct the arguments.
    ##########################################

    usage = "usage: %prog [options] <fasta input file> <fasta output file>\n\n" + \
            "    <fasta file>:  fasta sequence to annotate\n"
    description = "Ensure sequence names contain only alphanumeric characters\n"
    parser = OptionParser(usage=usage, description=description)

    options, args = parser.parse_args()

    if len(args) != 2:
        parser.print_help()
        return 1

    inputName = args[0]
    inputFile = open(inputName, "r")
    outputName = args[1]
    outputFile = open(outputName, "w")

    for header, seq in fastaRead(inputFile):
        fastaWrite(outputFile, fixHeader(header), seq)

    outputFile.close()
    inputFile.close()
    return 0

if __name__ == '__main__':
    exit(main())
