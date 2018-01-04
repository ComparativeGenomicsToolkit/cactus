#!/usr/bin/env python2.7
import os
import argparse



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--oldSamplingRates", type=str)
    parser.add_argument("--coveredSeedsFiles", type=str)
    parser.add_argument("--oldSamplingRate", type=float)
    parser.add_argument("--newSamplingRate", type=float)

    args = parser.parse_args()
    coveredSeedsFiles = args.coveredSeedsFiles.split(",")
    oldSamplingRates = args.oldSamplingRates
    
    seedCounts = {}
    for coveredSeedsFile in coveredSeedsFiles:
        with open(coveredSeedsFile, "r") as coveredSeedsFileRead:
            for line in coveredSeedsFileRead:
                lineinfo = line.split()
                if not len(lineinfo) == 2:
                    continue
                seed, count = lineinfo
                if seed in seedCounts:
                    seedCounts[seed] += count
                else:
                    seedCounts[seed] = count
    samplingRatesRead = open(oldSamplingRates, "r")
    for line in samplingRatesRead:
        seed, packed, samplingRate = line.split()
        samplingRate = float(samplingRate)
        if samplingRate < args.oldSamplingRate or seed in seedCounts:
            #Seed has been covered by an alignment in a previous iteration,
            #or in this one, so don't increment it
            print("%s\t%s\t%s" % (seed, packed, str(samplingRate)))
        else:
            print("%s\t%s\t%s" % (seed, packed, str(args.newSamplingRate)))
    samplingRatesRead.close()


if __name__ == "__main__":
    main()
