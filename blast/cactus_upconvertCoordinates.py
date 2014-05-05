#!/usr/bin/env python
from argparse import ArgumentParser
from collections import defaultdict
import sys
import os
from sonLib.bioio import cigarRead, cigarWrite, getTempFile, system

def getSequenceRanges(fa):
    """Get dict of (untrimmed header) -> [(range, trimmed header)] mappings
    from a trimmed fasta."""
    ret = defaultdict(list)
    curSeq = ""
    curHeader = None
    curTrimmedStart = None
    for line in fa:
        line = line.strip()
        if line == '':
            continue
        if line[0] == '>':
            if curHeader is not None:
                # Add previous seq info to dict
                trimmedRange = xrange(curTrimmedStart,
                                      curTrimmedStart + len(curSeq))
                untrimmedHeader = "|".join(curHeader.split("|")[:-1])
                ret[untrimmedHeader].append((trimmedRange, curHeader))
            curHeader = line[1:].split()[0]
            curTrimmedStart = int(curHeader.split('|')[-1])
            curSeq = ""
        else:
            curSeq += line
    if curHeader is not None:
        # Add final seq info to dict
        trimmedRange = xrange(curTrimmedStart,
                              curTrimmedStart + len(curSeq))
        untrimmedHeader = "|".join(curHeader.split("|")[:-1])
        ret[untrimmedHeader].append((trimmedRange, curHeader))
    for key in ret.keys():
        # Sort by range's start pos
        ret[key] = sorted(ret[key], key=lambda x: x[0][0])
    return ret

def validateRanges(seqRanges):
    """Fail if the given range dict contains overlapping ranges or if the
    ranges aren't sorted.
    """
    for seq, ranges in seqRanges.items():
        for i, range in enumerate(ranges):
            start = range[0][0]
            for j in xrange(i):
                range2 = ranges[j][0]
                assert start not in range2
                assert start >= range2[-1]
            for j in xrange(i+1,len(ranges)):
                range2 = ranges[j][0]
                assert start not in range2
                assert start < range2[0]

def sortCigarByContigAndPos(cigarPath, contigNum):
    contigNameKey = 2 if contigNum == 1 else 6
    startPosKey = 3 if contigNum == 1 else 7
    tempFile = getTempFile()
    system("sort -k %d -n %s | sort -k %d > %s" % (startPosKey, cigarPath, contigNameKey, tempFile))
    return tempFile

def upconvertCoords(cigarFile, seqRanges, contigNum):
    """Convert the coordinates of the given alignment, so that the
    alignment refers to a set of trimmed sequences originating from a
    contig rather than to the contig itself."""
    currentContig = None
    currentRangeIdx = None
    currentRange = None
    for alignment in cigarRead(cigarFile):
        contig = alignment.contig1 if contigNum == 1 else alignment.contig2
        if contig in seqRanges:
            minPos = min(alignment.start1, alignment.end1)
            maxPos = max(alignment.start1, alignment.end1)

            if contig != currentContig:
                currentContig = contig
                currentRangeIdx = 0
                currentRange = seqRanges[contig][0]
            while minPos not in currentRange[0] and currentRangeIdx < len(seqRanges[contig]) - 1:
                currentRangeIdx += 1
                currentRange = seqRanges[contig][currentRangeIdx]
            if minPos in currentRange[0]:
                if maxPos - 1 not in currentRange[0]:
                    raise RuntimeError("alignment on %s:%d-%d crosses "
                                       "trimmed sequence boundary" %\
                                       (alignment.contig1,
                                        alignment.start1,
                                        alignment.end1))
                if contigNum == 1:
                    alignment.start1 -= currentRange[0][0]
                    alignment.end1 -= currentRange[0][0]
                    alignment.contig1 = currentRange[1]
            else:
                raise RuntimeError("No trimmed sequence containing alignment "
                                   "on %s:%d-%d" % (alignment.contig1,
                                                    alignment.start1,
                                                    alignment.end1))
        cigarWrite(sys.stdout, alignment, False)

def main():
    parser = ArgumentParser()
    parser.add_argument("fasta", help="Trimmed fasta file")
    parser.add_argument("cigar", help="Alignments file to convert to trimmed "
                        "coordinates")
    parser.add_argument("contig", help="Contig # in to convert in each alignment (1 or 2)", type=int)
    args = parser.parse_args()
    assert args.contig == 1 or args.contig == 2
    seqRanges = getSequenceRanges(open(args.fasta))
    validateRanges(seqRanges)
    sortedCigarFile = sortCigarByContigAndPos(args.cigar, args.contig)
    upconvertCoords(open(sortedCigarFile), seqRanges, args.contig)
    os.remove(sortedCigarFile)

if __name__ == '__main__':
    main()
