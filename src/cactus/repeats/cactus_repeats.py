from cactus.shared.common import cactus_call, makeURL, RoundedJob
from cactus.shared.common import runGetChunks
from cactus.blast.blast import CollateBlasts

from toil.job import Job
from toil.common import Toil
from toil.lib.bioio import logger
from sonLib.bioio import getTempDirectory
from sonLib.bioio import catFiles

import argparse
import logging
import os


class CactusRepeatsOptions:
    def __init__(self):
        self.samplingRates = [0.01, 0.2, 1.0]
        self.chunkSize = 1000000

class GetInitialSamplingRates(RoundedJob):
    def __init__(self, sequenceIDs, samplingRate):
        RoundedJob.__init__(self)
        self.sequenceIDs = sequenceIDs
        self.samplingRate = samplingRate

    def run(self, fileStore):
        sequences = [fileStore.readGlobalFile(seqID) for seqID in self.sequenceIDs]
        combinedSeq = fileStore.getLocalTempFile()
        catFiles(sequences, combinedSeq)
        initialSamplingRates = fileStore.getLocalTempFile()
        with open(initialSamplingRates, "w") as fh:
            for line in cactus_call(check_output=True, parameters=["cPecanLastz",
                                "--tableonly=count",
                                "%s[unmask][multiple]" % os.path.basename(combinedSeq)]).split("\n"):
                lineinfo = line.split()
                if len(lineinfo) != 2:
                    continue
                seedInfo, count = line.split()
                packed, seed = seedInfo.split("/")
                #remove the colon
                seed = seed[:-1]
                fh.write("%s\t%s\t%f\n" % (seed, packed, float(self.samplingRate)))

        return fileStore.writeGlobalFile(initialSamplingRates)


class IterativeRepeatSampling(RoundedJob):
    def __init__(self, sequenceIDs, options):
        RoundedJob.__init__(self)
        self.sequenceIDs = sequenceIDs
        self.options = options

    def run(self, fileStore):
        sequences = [fileStore.readGlobalFile(seqID) for seqID in self.sequenceIDs]
        chunkDirectory = getTempDirectory(rootDir=fileStore.getLocalTempDir())
        chunks = runGetChunks(sequenceFiles=sequences, chunksDir=chunkDirectory,
                                   chunkSize=self.options.chunkSize,
                                   overlapSize=0)
        chunkIDs = [fileStore.writeGlobalFile(chunk) for chunk in chunks]

        initialSamplingRatesJob = GetInitialSamplingRates(sequenceIDs=self.sequenceIDs, samplingRate=self.options.samplingRates[0])

        prevJob = initialSamplingRatesJob
        samplingRatesID = self.addChild(initialSamplingRatesJob).rv()
        oldSamplingRate = self.options.samplingRates[0]
        for newSamplingRate in self.options.samplingRates[1:]:
            blastJob = RunBlastAndRaiseSamplingThreshold(chunkIDs=chunkIDs,
                    oldSamplingRate=oldSamplingRate, newSamplingRate=newSamplingRate,
                    options=self.options, samplingRatesID=samplingRatesID)
            prevJob.addFollowOn(blastJob)
            samplingRatesID = blastJob.rv()
            prevJob = blastJob
            oldSamplingRate = newSamplingRate
        return prevJob.rv()

class AdjustSamplingRates(RoundedJob):
    def __init__(self, coveredSeedsIDs, samplingRatesID, oldSamplingRate, newSamplingRate, options):
        RoundedJob.__init__(self)
        self.coveredSeedsIDs = coveredSeedsIDs
        self.samplingRatesID = samplingRatesID
        self.oldSamplingRate = oldSamplingRate
        self.newSamplingRate = newSamplingRate
        self.options = options
    def run(self, fileStore):
        logger.info("Raising sampling rate from %f to %f" % (self.oldSamplingRate, self.newSamplingRate))
        coveredSeedsFiles = [os.path.basename(fileStore.readGlobalFile(coveredSeedsID)) for coveredSeedsID in
                self.coveredSeedsIDs]
        samplingRates = fileStore.readGlobalFile(self.samplingRatesID)
        newSamplingRates = os.path.basename(fileStore.getLocalTempFile())
        logger.info("Work dir = %s" % os.path.dirname(samplingRates))
        cactus_call(outfile=newSamplingRates, parameters=["cactus_adjustSamplingRates",
                                                          "--oldSamplingRate", str(self.oldSamplingRate),
                                                          "--newSamplingRate", str(self.newSamplingRate),
                                                          "--oldSamplingRates", os.path.basename(samplingRates),
                                                          "--coveredSeedsFiles", ",".join(coveredSeedsFiles)])

        return fileStore.writeGlobalFile(newSamplingRates)


class RunBlastAndRaiseSamplingThreshold(RoundedJob):
    def __init__(self, chunkIDs, samplingRatesID, oldSamplingRate, newSamplingRate, options):
        RoundedJob.__init__(self)
        self.chunkIDs = chunkIDs
        self.samplingRatesID = samplingRatesID
        self.oldSamplingRate = oldSamplingRate
        self.newSamplingRate = newSamplingRate
        self.options = options
    def run(self, fileStore):
        return self.addChild(RunBlastAndRaiseSamplingThreshold2(chunkIDs=self.chunkIDs,
            samplingRatesID=self.samplingRatesID, oldSamplingRate=self.oldSamplingRate, newSamplingRate=self.newSamplingRate, options=self.options)).rv()

class RunBlastAndRaiseSamplingThreshold2(RoundedJob):
    def __init__(self, chunkIDs, samplingRatesID, oldSamplingRate, newSamplingRate, options):
        RoundedJob.__init__(self)
        self.chunkIDs = chunkIDs
        self.samplingRatesID = samplingRatesID
        self.oldSamplingRate = oldSamplingRate
        self.newSamplingRate = newSamplingRate
        self.options = options
    def run(self, fileStore):
        coveredSeedsIDs = []
        for i in range(len(self.chunkIDs)):
            for j in range(i, len(self.chunkIDs)):
                coveredSeedsIDs.append(self.addChild(RunBlast(options=self.options, chunkID1=self.chunkIDs[i], chunkID2=self.chunkIDs[j], samplingRatesID=self.samplingRatesID)).rv())
        return self.addFollowOn(AdjustSamplingRates(coveredSeedsIDs=coveredSeedsIDs,
            samplingRatesID=self.samplingRatesID, oldSamplingRate=self.oldSamplingRate, newSamplingRate=self.newSamplingRate, options=self.options)).rv()


class RunBlast(RoundedJob):
    def __init__(self, options, chunkID1, chunkID2, samplingRatesID):
        RoundedJob.__init__(self)
        self.options = options
        self.chunkID1 = chunkID1
        self.chunkID2 = chunkID2
        self.samplingRatesID = samplingRatesID
    def run(self, fileStore):
        samplingRates = fileStore.readGlobalFile(self.samplingRatesID)
        chunk1 = fileStore.readGlobalFile(self.chunkID1)
        chunk2 = fileStore.readGlobalFile(self.chunkID2)
        lastZSequenceHandling  = ['%s[multiple,unmask][nameparse=darkspace]' % os.path.basename(chunk1), '%s[nameparse=darkspace][unmask]' % os.path.basename(chunk2)]
        alignments = fileStore.getLocalTempFile()
        logger.info("Work dir = %s" % os.path.dirname(chunk1))
        cactus_call(outfile=alignments,
                    parameters=["cPecanLastz"] + lastZSequenceHandling +
                                ["--samplingRates=%s" % os.path.basename(samplingRates),
                                "--notrivial",
                                "--format=general:name1,zstart1,end1,name2,zstart2+,end2+", "--markend"])

        coveredSeeds = fileStore.getLocalTempFile()
        cactus_call(outfile=coveredSeeds,
                    parameters=["cactus_coveredSeeds",
                                "--seq", chunk1,
                                "--alignments", alignments])
        return fileStore.writeGlobalFile(coveredSeeds)

def main():
    parser = argparse.ArgumentParser()
    Job.Runner.addToilOptions(parser)
    parser.add_argument("--sequences")
    parser.add_argument("--samplingRates")
    parser.add_argument("--outfile")

    args = parser.parse_args()

    args.disableCaching = "True"

    with Toil(args) as toil:
        sequenceIDs = [toil.importFile(makeURL(seq)) for seq in args.sequences.split(",")]
        repeatsOptions = CactusRepeatsOptions()
        repeatsOptions.samplingRates = args.samplingRates.split(",")
        repeatsOptions.samplingRates = [float(rate) for rate in repeatsOptions.samplingRates]
        samplingRatesID = toil.start(IterativeRepeatSampling(sequenceIDs, repeatsOptions))
        toil.exportFile(samplingRatesID, makeURL(args.outfile))

if __name__ == "__main__":
    main()
