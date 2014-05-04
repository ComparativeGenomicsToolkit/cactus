#!/usr/bin/env python
from argparse import ArgumentParser
from collections import defaultdict
import sys
from sonLib.bioio import cigarRead, cigarWrite

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
    return ret

def validateRanges(seqRanges):
    """Fail if the given range dict contains overlapping ranges."""
    for seq, ranges in seqRanges.items():
        for i, range in enumerate(ranges):
            start = range[0]
            for j in xrange(i):
                range2 = ranges[j]
                assert start not in range2
            for j in xrange(i+1,len(ranges)):
                range2 = ranges[j]
                assert start not in range2

def upconvertCoords(cigarFile, seqRanges):
    for alignment in cigarRead(cigarFile):
        if alignment.contig1 not in seqRanges and alignment.contig2 not in seqRanges:
            # alignment is not relevant to the trimmed fasta, leave it as is
            cigarWrite(stdout, alignment, False)
            continue
        if alignment.contig1 in seqRanges:
            # Search through available ranges--could be a binary
            # search if this is slow.
            found = False
            for range in seqRanges[alignment.contig1]:
                minPos = min(alignment.start1, alignment.end1)
                maxPos = max(alignment.start1, alignment.end1)
                if minPos in range[0]:
                    if maxPos - 1 not in range[0]:
                        raise RuntimeError("alignment on %s:%d-%d crosses "
                                           "trimmed sequence boundary" %\
                                           (alignment.contig1,
                                            alignment.start1,
                                            alignment.end1))
                    alignment.start1 -= range[0][0]
                    alignment.end1 -= range[0][0]
                    alignment.contig1 = alignment.contig1 + "|%d" % range[0][0]
                    found = True
                    break
            if not found:
                raise RuntimeError("No trimmed sequence containing alignment "
                                   "on %s:%d-%d" % (alignment.contig1,
                                                    alignment.start1,
                                                    alignment.end1))
        if alignment.contig2 in seqRanges:
            # Search through available ranges--could be a binary
            # search if this is slow.
            found = False
            for range in seqRanges[alignment.contig2]:
                minPos = min(alignment.start2, alignment.end2)
                maxPos = max(alignment.start2, alignment.end2)
                if minPos in range[0]:
                    if maxPos - 1 not in range[0]:
                        raise RuntimeError("alignment on %s:%d-%d crosses "
                                           "trimmed sequence boundary" %\
                                           (alignment.contig2,
                                            alignment.start2,
                                            alignment.end2))
                    alignment.start2 -= range[0][0]
                    alignment.end2 -= range[0][0]
                    alignment.contig2 = alignment.contig2 + "|%d" % range[0][0]
                    found = True
                    break
            if not found:
                raise RuntimeError("No trimmed sequence containing alignment "
                                   "on %s:%d-%d" % (alignment.contig2,
                                                    alignment.start2,
                                                    alignment.end2))
        cigarWrite(sys.stdout, alignment, False)

def main():
    parser = ArgumentParser()
    parser.add_argument("fasta", help="Trimmed fasta file")
    parser.add_argument("cigar", help="Alignments file to convert to trimmed "
                        "coordinates")
    args = parser.parse_args()
    seqRanges = getSequenceRanges(open(args.fasta))
    validateRanges(seqRanges)
    upconvertCoords(open(args.cigar), seqRanges)

if __name__ == '__main__':
    main()
