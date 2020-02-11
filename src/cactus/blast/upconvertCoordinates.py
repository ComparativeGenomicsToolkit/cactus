#!/usr/bin/env python3
from argparse import ArgumentParser
from collections import defaultdict
import sys
import os
from sonLib.bioio import cigarRead, cigarWrite, getTempFile, system

def getSequenceRanges(fa):
    """Get dict of (untrimmed header) -> [(start, non-inclusive end)] mappings
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
                trimmedRange = (curTrimmedStart,
                                curTrimmedStart + len(curSeq))
                untrimmedHeader = "|".join(curHeader.split("|")[:-1])
                ret[untrimmedHeader].append(trimmedRange)
            curHeader = line[1:].split()[0]
            curTrimmedStart = int(curHeader.split('|')[-1])
            curSeq = ""
        else:
            curSeq += line
    if curHeader is not None:
        # Add final seq info to dict
        trimmedRange = (curTrimmedStart,
                        curTrimmedStart + len(curSeq))
        untrimmedHeader = "|".join(curHeader.split("|")[:-1])
        ret[untrimmedHeader].append(trimmedRange)
    for key in list(ret.keys()):
        # Sort by range's start pos
        ret[key] = sorted(ret[key], key=lambda x: x[0])
    return ret

def validateRanges(seqRanges):
    """Fail if the given range dict contains overlapping ranges or if the
    ranges aren't sorted.
    """
    for seq, ranges in list(seqRanges.items()):
        for i, range in enumerate(ranges):
            start = range[0]
            if i - 1 >= 0:
                range2 = ranges[i - 1]
                assert start >= range2[1]
            if i + 1 < len(ranges):
                range2 = ranges[i + 1]
                assert start < range2[0]

def sortCigarByContigAndPos(cigarPath, contigNum):
    contigNameKey = 2 if contigNum == 1 else 6
    startPosKey = 3 if contigNum == 1 else 7
    tempFile = getTempFile()
    system("sort -k %d,%d -k %d,%dn %s > %s" % (contigNameKey, contigNameKey, startPosKey, startPosKey, cigarPath, tempFile))
    return tempFile

def upconvertCoords(cigarPath, fastaPath, contigNum, outputFile):
    """Convert the coordinates of the given alignment, so that the
    alignment refers to a set of trimmed sequences originating from a
    contig rather than to the contig itself."""
    with open(fastaPath) as f:
        seqRanges = getSequenceRanges(f)
    validateRanges(seqRanges)
    sortedCigarPath = sortCigarByContigAndPos(cigarPath, contigNum)
    sortedCigarFile = open(sortedCigarPath)

    currentContig = None
    currentRangeIdx = None
    currentRange = None
    for alignment in cigarRead(sortedCigarFile):
        # contig1 and contig2 are reversed in python api!!
        contig = alignment.contig2 if contigNum == 1 else alignment.contig1
        minPos = min(alignment.start2, alignment.end2) if contigNum == 1 else min(alignment.start1, alignment.end1)
        maxPos = max(alignment.start2, alignment.end2) if contigNum == 1 else max(alignment.start1, alignment.end1)
        if contig in seqRanges:
            if contig != currentContig:
                currentContig = contig
                currentRangeIdx = 0
                currentRange = seqRanges[contig][0]
            while (minPos >= currentRange[1] or minPos < currentRange[0]) and currentRangeIdx < len(seqRanges[contig]) - 1:
                currentRangeIdx += 1
                currentRange = seqRanges[contig][currentRangeIdx]
            if currentRange[0] <= minPos < currentRange[1]:
                if maxPos - 1 > currentRange[1]:
                    raise RuntimeError("alignment on %s:%d-%d crosses "
                                       "trimmed sequence boundary" %\
                                       (contig,
                                        minPos,
                                        maxPos))
                if contigNum == 1:
                    alignment.start2 -= currentRange[0]
                    alignment.end2 -= currentRange[0]
                    alignment.contig2 = contig + ("|%d" % currentRange[0])
                else:
                    alignment.start1 -= currentRange[0]
                    alignment.end1 -= currentRange[0]
                    alignment.contig1 = contig + ("|%d" % currentRange[0])
            else:
                raise RuntimeError("No trimmed sequence containing alignment "
                                   "on %s:%d-%d" % (contig,
                                                    minPos,
                                                    maxPos))
        cigarWrite(outputFile, alignment, False)
    os.remove(sortedCigarPath)
