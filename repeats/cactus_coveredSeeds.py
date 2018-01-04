#!/usr/bin/env python

import argparse

from sonLib.bioio import fastaRead

def getSeed(seedPattern, kmer):
    seed = ""
    for base, keep in zip(kmer, seedPattern):
        if bool(int(keep)):
            seed += base.upper()
        else:
            seed += 'x'
    return seed

def iterSeeds(seq, start, end, seedPattern):
    kmerLength = len(seedPattern)
    for i in range(start, end - kmerLength):
        yield getSeed(seedPattern=seedPattern, kmer=seq[i: i + kmerLength])
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seq", type=str)
    parser.add_argument("--alignments", type=str)
    parser.add_argument("--seedPattern", type=str, default="1110100110010101111")

    args = parser.parse_args()

    
    alignmentsFile = open(args.alignments,"rt")

    chromToIntervals = {}

    lineNumber = 0
    for line in alignmentsFile:
        lineNumber += 1
        line = line.strip()
        if (line == "") or (line.startswith("#")): continue

        fields = line.split()
        assert (len(fields) >= 3), \
              "not enough fields (line %s): %s" % (lineNumber,line)

        try:
                chrom  = fields[0]
                start = int(fields[1])
                end   = int(fields[2])
                if (start < 0):    raise ValueError
                if (start >= end): raise ValueError
        except ValueError:
                assert (False), \
                      "bad line (line %s): %s" % (lineNumber,line)

        if (chrom not in chromToIntervals): chromToIntervals[chrom] = []
        chromToIntervals[chrom] += [(start,end)]

    alignmentsFile.close()
    print("Found chromosomes %s" % chromToIntervals.keys())

    seqFile = open(args.seq, "r")
    seenSeeds = {}

    for (chrom, seq) in fastaRead(seqFile):
        if not chrom in chromToIntervals:
            continue
        for (start, end) in chromToIntervals[chrom]:
            for seed in iterSeeds(start=start, end=end, seq=seq, seedPattern=args.seedPattern):
                if seed in seenSeeds:
                    seenSeeds[seed] += 1
                else:
                    seenSeeds[seed] = 1
    seqFile.close()
    
    for seed in seenSeeds:
        print "%s %i" % (seed, seenSeeds[seed])

if __name__ == "__main__":
    main()
