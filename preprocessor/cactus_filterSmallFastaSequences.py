#!/usr/bin/env python3

## Filter small sequences out of a fasta file.  For use with flies,
## for example, where scaffolds of length <200kb seem to be considered
## no mans land

import os
from optparse import OptionParser

from sonLib.bioio import fastaRead
from sonLib.bioio import fastaWrite
from sonLib.bioio import getTempFile


# for every sequence, determine if its contained in the file
# (starts with |1|0; and there is a differently named sequence after it),
# and its length (defined as the max offset + length) found for the name
# assumption: sequences with same name are contiguous (which is true for
# cactus_batchChunk output, which this is tailored for)
# **only bother if names seem to be in chunk format (None returned otherwise)
def containedSequences(inputFile):
    lookup = dict()
    prev = ""
    for header, seq in fastaRead(inputFile):
        if '|1|' not in header:
            assert len(lookup) == 0
            return None
        else:
            idx = header.find('|1|')
            name = header[:idx]
            offset = header[idx+3:]
            if offset.isdigit() == False:
                assert len(lookup) == 0
                return None
            if int(offset) == 0:
                assert (name in lookup) == False
                lookup[name] = (len(seq), False)
            elif (name in lookup) == True:
                lookup[name] = (max(lookup[name][0], int(offset) + len(seq)), lookup[name][1])
            if name != prev and prev in lookup:
                lookup[prev] = (lookup[prev][0], True)
            prev = name
    return lookup

def tooShort(header, seq, options, contTable):
    isTooShort = False
    if contTable is not None:
        key = header[:header.find('|1|')]
        if key in contTable:
            length, flag = contTable[key]
            isTooShort = flag and length < options.length
    else:
        isTooShort = len(seq) < options.length

    return isTooShort

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

    contTable = containedSequences(inputFile)
    inputFile.seek(0)

    for header, seq in fastaRead(inputFile):
        if tooShort(header, seq, options, contTable) == False:
            fastaWrite(outputFile, header, seq)

    outputFile.close()
    inputFile.close()
    return 0

if __name__ == '__main__':
    exit(main())
