#!/usr/bin/env python3

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import os
import sys
import random
import time
import shutil
import filecmp
import pytest

from sonLib.bioio import system
from sonLib.bioio import logger
from sonLib.bioio import fastaWrite
from sonLib.bioio import getRandomSequence
from sonLib.bioio import mutateSequence
from sonLib.bioio import reverseComplement
from sonLib.bioio import getTempFile
from sonLib.bioio import getTempDirectory
from sonLib.bioio import cigarRead
from sonLib.bioio import PairwiseAlignment
from sonLib.bioio import TestStatus
from sonLib.bioio import catFiles
from sonLib.bioio import popenCatch

from cactus.shared.test import checkCigar, getTestLogLevel
from cactus.blast.blast import decompressFastaFile, compressFastaFile

from cactus.shared.common import runLastz
from cactus.shared.common import makeURL

from cactus.blast.blast import BlastOptions
from cactus.blast.blast import BlastIngroupsAndOutgroups
from cactus.blast.blast import BlastSequencesAllAgainstAll
from cactus.blast.blast import BlastSequencesAgainstEachOther
from cactus.blast.blast import calculateCoverage

from toil.job import Job
from toil.common import Toil

@pytest.mark.blast
@TestStatus.needsTestData
class TestCase(unittest.TestCase):
    def setUp(self):
        self.testNo = TestStatus.getTestSetup(1, 5, 10, 100)
        self.tempDir = getTempDirectory(os.getcwd())
        self.tempFiles = []
        unittest.TestCase.setUp(self)
        self.tempOutputFile = os.path.join(self.tempDir, "results1.txt")
        self.tempFiles.append(self.tempOutputFile)
        self.tempOutputFile2 = os.path.join(self.tempDir, "results2.txt")
        self.tempFiles.append(self.tempOutputFile2)
        self.encodePath = os.path.join(TestStatus.getPathToDataSets(), "MAY-2005")

    def tearDown(self):
        for tempFile in self.tempFiles:
            if os.path.exists(tempFile):
                os.remove(tempFile)
        unittest.TestCase.tearDown(self)
        system("rm -rf %s" % self.tempDir)

    def runComparisonOfBlastScriptVsNaiveBlast(self, blastMode):
        """We compare the output with a naive run of the blast program, to check the results are nearly
        equivalent.
        """
        encodeRegions = [ "ENm00" + str(i) for i in range(1,2) ] #, 2) ] #Could go to six
        species = ("human", "mouse", "dog")
        #Other species to try "rat", "monodelphis", "macaque", "chimp"
        for encodeRegion in encodeRegions:
            regionPath = os.path.join(self.encodePath, encodeRegion)
            for i in range(len(species)):
                species1 = species[i]
                for species2 in species[i+1:]:
                    seqFile1 = os.path.join(regionPath, "%s.%s.fa" % (species1, encodeRegion))
                    seqFile2 = os.path.join(regionPath, "%s.%s.fa" % (species2, encodeRegion))

                    #Run simple blast
                    runNaiveBlast(seqFile1, seqFile2, self.tempOutputFile, self.tempDir)
                    logger.info("Ran the naive blast okay")

                    #Run cactus blast pipeline
                    toilDir = os.path.join(getTempDirectory(self.tempDir), "toil")
                    if blastMode == "allAgainstAll":
                        runCactusBlast(sequenceFiles=[ seqFile1, seqFile2 ],
                                       alignmentsFile=self.tempOutputFile2, toilDir=toilDir,
                                       chunkSize=500000, overlapSize=10000)
                    else:
                        runCactusBlast(sequenceFiles=[ seqFile1 ], alignmentsFile=self.tempOutputFile2,
                                       toilDir=toilDir, chunkSize=500000, overlapSize=10000,
                                       targetSequenceFiles=[ seqFile2 ])
                    logger.info("Ran cactus_blast okay")
                    logger.critical("Comparing cactus_blast and naive blast; using mode: %s" % blastMode)
                    checkCigar(self.tempOutputFile)
                    checkCigar(self.tempOutputFile2)
                    compareResultsFile(self.tempOutputFile, self.tempOutputFile2)

    @TestStatus.mediumLength
    def testBlastEncodeAllAgainstAll(self):
        """For each encode region, for set of pairwise species, run
        cactus_blast.py in all-against-all mode.
        """
        self.runComparisonOfBlastScriptVsNaiveBlast(blastMode="allAgainstAll")

    @TestStatus.mediumLength
    def testBlastEncode(self):
        """For each encode region, for set of pairwise species, run
        cactus_blast.py in one set of sequences against another set mode.
        """
        self.runComparisonOfBlastScriptVsNaiveBlast(blastMode="againstEachOther")

    @TestStatus.longLength
    def testAddingOutgroupsImprovesResult(self):
        """Run blast on "ingroup" and "outgroup" encode regions, and ensure
        that adding an extra outgroup only adds alignments if
        possible, and doesn't lose any
        """
        encodeRegion = "ENm001"
        ingroups = ["human", "macaque"]
        outgroups = ["rabbit", "dog", "rat", "platypus", "xenopus", "fugu"]
        MAX_NUM_OUTGROUPS = 3
        # subselect a random set of outgroups in the same order
        outgroups = [outgroups[i] for i in sorted(random.sample(list(range(len(outgroups))), MAX_NUM_OUTGROUPS))]
        regionPath = os.path.join(self.encodePath, encodeRegion)
        ingroupPaths = [os.path.join(regionPath, x + "." + encodeRegion + ".fa") for x in ingroups]
        outgroupPaths = [os.path.join(regionPath, x + "." + encodeRegion + ".fa") for x in outgroups]
        results = []
        for numOutgroups in range(1, len(outgroups) + 1):
            # Align w/ increasing numbers of outgroups
            subResults = getTempFile()
            subOutgroupPaths = outgroupPaths[:numOutgroups]
            print(("aligning %s vs %s" % (",".join(ingroupPaths), ",".join(subOutgroupPaths))))
            tmpToil = os.path.join(self.tempDir, "outgroupToil")
            runCactusBlastIngroupsAndOutgroups(ingroupPaths, subOutgroupPaths, alignmentsFile=subResults, toilDir=tmpToil)
            results.append(subResults)

        # Print diagnostics about coverage
        for i, subResults in enumerate(results):
            for ingroup, ingroupPath in zip(ingroups, ingroupPaths):
                ingroupCoverage = getTempFile(rootDir=self.tempDir)
                calculateCoverage(sequenceFile=ingroupPath, cigarFile=subResults, outputFile=ingroupCoverage)
                coveredBases = popenCatch("cat %s | awk '{ total += $3 - $2 } END { print total }'" % ingroupCoverage)
                print(("covered bases on %s using %d outgroups: %s" % (ingroup, i + 1, coveredBases)))

        resultsSets = [loadResults(x) for x in results]
        for i, moreOutgroupsResults in enumerate(resultsSets[1:]):
            # Make sure the results from (n+1) outgroups are
            # (very nearly) a superset of the results from n outgroups
            print(("Using %d addl outgroup(s):" % (i + 1)))
            comparator =  ResultComparator(resultsSets[0], moreOutgroupsResults)
            print(comparator)
            self.assertTrue(comparator.sensitivity >= 0.99)

        # Ensure that the new alignments don't cover more than
        # x% of already existing alignments to human
        for i in range(1, len(resultsSets)):
            prevResults = resultsSets[i-1][0]
            curResults = resultsSets[i][0]
            prevResultsHumanPos = set([(x[0], x[1]) if "human" in x[0] else (x[2], x[3]) for x in [x for x in prevResults if "human" in x[0] or "human" in x[2]]])
            newAlignments = curResults.difference(prevResults)
            newAlignmentsHumanPos =  set([(x[0], x[1]) if "human" in x[0] else (x[2], x[3]) for x in [x for x in newAlignments if "human" in x[0] or "human" in x[2]]])
            print(("addl outgroup %d:" % i))
            print(("bases re-covered: %f (%d)" % (len(newAlignmentsHumanPos.intersection(prevResultsHumanPos))/float(len(prevResultsHumanPos)), len(newAlignmentsHumanPos.intersection(prevResultsHumanPos)))))
        for subResult in results:
            os.remove(subResult)

    @TestStatus.mediumLength
    def testKeepingCoverageOnIngroups(self):
        """Tests whether the --ingroupCoverageDir option works as
        advertised."""
        encodeRegion = "ENm001"
        ingroups = ["human", "cow"]
        outgroups = ["macaque", "rabbit", "dog"]
        regionPath = os.path.join(self.encodePath, encodeRegion)
        ingroupPaths = [os.path.join(regionPath, x + "." + encodeRegion + ".fa") for x in ingroups]
        outgroupPaths = [os.path.join(regionPath, x + "." + encodeRegion + ".fa") for x in outgroups]
        # Run blast in "ingroup vs outgroups" mode, requesting to keep
        # the bed files that show outgroup coverage on the ingroup.
        toilDir = os.path.join(self.tempDir, "tmp_toil")
        outgroupFragmentPaths = [getTempFile(rootDir=self.tempDir) for outgroup in outgroups]
        ingroupCoveragePaths = [getTempFile(rootDir=self.tempDir) for ingroup in ingroups]
        runCactusBlastIngroupsAndOutgroups(ingroups=ingroupPaths, outgroups=outgroupPaths, alignmentsFile=self.tempOutputFile, outgroupFragmentPaths=outgroupFragmentPaths, ingroupCoveragePaths=ingroupCoveragePaths, toilDir=toilDir)
        for i, ingroupPath in enumerate(ingroupPaths):
            # Get the coverage from the outgroups independently and
            # check that it's the same as the file in
            # ingroupCoverageDir
            otherIngroupPath = ingroupPaths[1] if i == 0 else ingroupPaths[0]
            # To filter out alignments from the other ingroup and
            # self-alignments we need to create a fasta with all the
            # outgroup fragments in it.
            outgroupsCombined = getTempFile(rootDir=self.tempDir)
            for outgroupFragmentPath in outgroupFragmentPaths:
                system("cat %s >> %s" % (outgroupFragmentPath, outgroupsCombined))
            independentCoverageFile = getTempFile(rootDir=self.tempDir)
            calculateCoverage(fromGenome=outgroupsCombined, sequenceFile=ingroupPath, cigarFile=self.tempOutputFile, outputFile=independentCoverageFile)

            # find the coverage file cactus_blast kept (should be
            # named according to the basename of the ingroup path
            # file)
            keptCoverageFile = ingroupCoveragePaths[i]
            self.assertTrue(filecmp.cmp(independentCoverageFile, keptCoverageFile))

    @TestStatus.mediumLength
    def testProgressiveOutgroupsVsAllOutgroups(self):
        """Tests the difference in outgroup coverage on an ingroup when
        running in "ingroups vs. outgroups" mode and "set against set"
        mode.
        """
        encodeRegion = "ENm001"
        ingroup = "human"
        outgroups = ["macaque", "rabbit", "dog"]
        regionPath = os.path.join(self.encodePath, encodeRegion)
        ingroupPath = os.path.join(regionPath, ingroup + "." + encodeRegion + ".fa")
        outgroupPaths = [os.path.join(regionPath, x + "." + encodeRegion + ".fa") for x in outgroups]
        # Run in "set against set" mode, aligning the entire ingroup
        # vs each outgroup
        runCactusBlast([ingroupPath], alignmentsFile=self.tempOutputFile,
                       toilDir=os.path.join(self.tempDir, "setVsSetToil"),
                       chunkSize=500000, overlapSize=10000,
                       targetSequenceFiles=outgroupPaths)
        # Run in "ingroup vs outgroups" mode, aligning the ingroup vs
        # the outgroups in order, trimming away sequence that's
        # already been aligned.
        runCactusBlastIngroupsAndOutgroups([ingroupPath], outgroupPaths, alignmentsFile=self.tempOutputFile2, toilDir=os.path.join(self.tempDir, "outgroupToil"))

        # Get the coverage on the ingroup, in bases, from each run.
        coverageSetVsSetUnfiltered = getTempFile(rootDir=self.tempDir)
        calculateCoverage(sequenceFile=ingroupPath, cigarFile=self.tempOutputFile, outputFile=coverageSetVsSetUnfiltered)
        coverageSetVsSet = int(popenCatch("cat %s | awk '{ total +=  $3 - $2} END { print total }'" % coverageSetVsSetUnfiltered))
        coverageIngroupVsOutgroupsUnfiltered = getTempFile(rootDir=self.tempDir)
        calculateCoverage(sequenceFile=ingroupPath, cigarFile=self.tempOutputFile2, outputFile=coverageIngroupVsOutgroupsUnfiltered)
        coverageIngroupVsOutgroups = int(popenCatch("cat %s | awk '{ total +=  $3 - $2} END { print total }'" % coverageIngroupVsOutgroupsUnfiltered))

        print(("total coverage on human (set vs set mode, %d outgroups): %d" % (len(outgroups), coverageSetVsSet)))
        print(("total coverage on human (ingroup vs outgroup mode, %d outgroups): %d" % (len(outgroups), coverageIngroupVsOutgroups)))

        # Make sure we're getting a reasonable fraction of the
        # alignments when using the trimming strategy.
        self.assertTrue(float(coverageIngroupVsOutgroups)/coverageSetVsSet >= 0.95)

        # Get the coverage on the ingroup, in bases, from just the
        # last outgroup. Obviously this should be much higher in set
        # vs set mode than in ingroup vs outgroup mode.
        outgroupAlignments = getTempFile(rootDir=self.tempDir)
        system("grep %s %s > %s" % (outgroups[-1], self.tempOutputFile, outgroupAlignments))
        coverageFileSetVsSet = getTempFile(rootDir=self.tempDir)
        calculateCoverage(sequenceFile=ingroupPath, cigarFile=outgroupAlignments, outputFile=coverageFileSetVsSet)

        coverageFromLastOutgroupSetVsSet = int(popenCatch("cat %s | awk '{ total +=  $3 - $2} END { print total }'" % coverageFileSetVsSet))


        outgroupAlignments = getTempFile(rootDir=self.tempDir)
        system("grep %s %s > %s" % (outgroups[-1], self.tempOutputFile2, outgroupAlignments))
        coverageFileInVsOut = getTempFile(rootDir=self.tempDir)
        calculateCoverage(sequenceFile=ingroupPath, cigarFile=outgroupAlignments, outputFile=coverageFileInVsOut)
        coverageFromLastOutgroupInVsOut = int(popenCatch("cat %s | awk '{ total +=  $3 - $2} END { print total }'" % coverageFileInVsOut))

        print(("total coverage on human from last outgroup in set (%s) (set vs set mode): %d" % (outgroups[-1], coverageFromLastOutgroupSetVsSet)))
        print(("total coverage on human from last outgroup in set (%s) (ingroup vs outgroup mode): %d" % (outgroups[-1], coverageFromLastOutgroupInVsOut)))

        self.assertTrue(float(coverageFromLastOutgroupInVsOut)/coverageFromLastOutgroupSetVsSet <= 0.10)

    @TestStatus.mediumLength
    def testBlastParameters(self):
        """Tests if changing parameters of lastz creates results similar to the desired default.
        """
        encodeRegion = "ENm001"
        species = ("human", "mouse", "dog")
        #Other species to try "rat", "monodelphis", "macaque", "chimp"
        regionPath = os.path.join(self.encodePath, encodeRegion)
        for i in range(len(species)):
            species1 = species[i]
            for species2 in species[i+1:]:
                seqFile1 = os.path.join(regionPath, "%s.%s.fa" % (species1, encodeRegion))
                seqFile2 = os.path.join(regionPath, "%s.%s.fa" % (species2, encodeRegion))

                #Run the random
                runNaiveBlast(seqFile1, seqFile2, self.tempOutputFile2, self.tempDir,
                              lastzArguments="--nogapped --hspthresh=3000 --ambiguous=iupac")
                #Run the blast
                runNaiveBlast(seqFile1, seqFile2, self.tempOutputFile, self.tempDir,
                              lastzArguments="--nogapped --step=3 --hspthresh=3000 --ambiguous=iupac")

                logger.critical("Comparing blast settings")
                compareResultsFile(self.tempOutputFile, self.tempOutputFile2, 0.7)

    @unittest.skip("fails on no chunks produced error")
    @TestStatus.mediumLength
    def testBlastRandom(self):
        """Make some sequences, put them in a file, call blast with random parameters
        and check it runs okay.
        """
        tempSeqFile = os.path.join(self.tempDir, "tempSeq.fa")
        self.tempFiles.append(tempSeqFile)
        for test in range(self.testNo):
            seqNo = random.choice(list(range(0, 10)))
            seq = getRandomSequence(8000)[1]
            fileHandle = open(tempSeqFile, 'w')
            for fastaHeader, seq in [ (str(i), mutateSequence(seq, 0.3*random.random())) for i in range(seqNo) ]:
                if random.random() > 0.5:
                    seq = reverseComplement(seq)
                fastaWrite(fileHandle, fastaHeader, seq)
            fileHandle.close()
            chunkSize = random.choice(list(range(500, 9000)))
            overlapSize = random.choice(list(range(2, 100)))
            toilDir = os.path.join(getTempDirectory(self.tempDir), "toil")
            runCactusBlast([ tempSeqFile ], self.tempOutputFile, toilDir, chunkSize, overlapSize)
            system("rm -rf %s " % toilDir)

    @TestStatus.mediumLength
    def testCompression(self):
        tempSeqFile = os.path.join(self.tempDir, "tempSeq.fa")
        tempSeqFile2 = os.path.join(self.tempDir, "tempSeq2.fa")
        self.tempFiles.append(tempSeqFile)
        self.tempFiles.append(tempSeqFile2)
        self.encodePath = os.path.join(self.encodePath, "ENm001")
        catFiles([ os.path.join(self.encodePath, fileName) for fileName in os.listdir(self.encodePath) ], tempSeqFile)
        startTime = time.time()
        compressFastaFile(tempSeqFile)
        logger.critical("It took %s seconds to compress the fasta file" % (time.time() - startTime))
        startTime = time.time()
        system("rm %s" % tempSeqFile + ".bz2")
        system("bzip2 --keep --fast %s" % tempSeqFile)
        logger.critical("It took %s seconds to compress the fasta file by system functions" % (time.time() - startTime))
        startTime = time.time()
        decompressFastaFile(tempSeqFile + ".bz2", tempSeqFile2)
        logger.critical("It took %s seconds to decompress the fasta file" % (time.time() - startTime))
        system("rm %s" % tempSeqFile2)
        startTime = time.time()
        system("bunzip2 --stdout %s > %s" % (tempSeqFile + ".bz2", tempSeqFile2))
        logger.critical("It took %s seconds to decompress the fasta file using system function" % (time.time() - startTime))
        logger.critical("File sizes, before: %s, compressed: %s, decompressed: %s" % (os.stat(tempSeqFile).st_size, os.stat(tempSeqFile + ".bz2").st_size, os.stat(tempSeqFile2).st_size))
        #Above test justifies out use of compression to reduce network transfer!
        #startTime = time.time()
        #runNaiveBlast([ tempSeqFile ], self.tempOutputFile, self.tempDir, lastzOptions="--nogapped --step=3 --hspthresh=3000 --ambiguous=iupac")
        #logger.critical("It took %s seconds to run blast" % (time.time() - startTime))


def compareResultsFile(results1, results2, closeness=0.95):
    results1 = loadResults(results1)
    logger.info("Loaded first results")

    #Now compare the results
    results2 = loadResults(results2)
    logger.info("Loaded second results")

    checkResultsAreApproximatelyEqual(ResultComparator(results2, results1), closeness)
    logger.info("Compared the blast results, using the first results as the 'true' results, and the second results as the predicted results")


def checkResultsAreApproximatelyEqual(resultsComparator, closeness=0.95):
    """Checks that the comparisons show the two blasts are suitably close together.
    """
    logger.critical("Results are: %s" % resultsComparator)
    assert resultsComparator.sensitivity >= closeness
    assert resultsComparator.specificity >= closeness

class ResultComparator:
    def __init__(self, trueResults, predictedResults):
        """Compares two sets of results and returns a set of statistics comparing them.
        """
        #Totals
        self.trueHits = trueResults[1]
        self.trueLength = len(trueResults[0])
        self.predictedHits = predictedResults[1]
        self.predictedLength = len(predictedResults[0])
        self.intersectionSize = len(trueResults[0].intersection(predictedResults[0]))
        #Sensitivity
        self.trueDifference = float(self.trueLength - self.intersectionSize)
        self.sensitivity = 1.0 - self.trueDifference / self.trueLength
        #Specificity
        self.predictedDifference = float(self.predictedLength - self.intersectionSize)
        self.specificity = 1.0 - self.predictedDifference / self.predictedLength
        #Union size
        self.unionSize = self.intersectionSize + self.trueDifference + self.predictedDifference
        #Symmetric difference
        self.symmDiff = self.intersectionSize / self.unionSize

    def __str__(self):
        return "True length: %s, predicted length: %s, union size: %s, \
intersection size: %s, symmetric difference: %s, \
true difference: %s, sensitivity: %s, \
predicted difference: %s, specificity: %s" % \
    (self.trueLength, self.predictedLength, self.unionSize,
     self.intersectionSize, self.symmDiff,
     self.trueDifference, self.sensitivity,
     self.predictedDifference, self.specificity)

def loadResults(resultsFile):
    """Puts the results in a set.
    """
    pairsSet = set()
    fileHandle = open(resultsFile, 'r')
    totalHits = 0
    for pairwiseAlignment in cigarRead(fileHandle):
        totalHits +=1
        i = pairwiseAlignment.start1
        s1 = 1
        if not pairwiseAlignment.strand1:
            i -= 1
            s1 = -1

        j = pairwiseAlignment.start2
        s2 = 1
        if not pairwiseAlignment.strand2:
            j -= 1
            s2 = -1

        for operation in pairwiseAlignment.operationList:
            if operation.type == PairwiseAlignment.PAIRWISE_INDEL_X:
                i += operation.length * s1
            elif operation.type == PairwiseAlignment.PAIRWISE_INDEL_Y:
                j += operation.length * s2
            else:
                assert operation.type == PairwiseAlignment.PAIRWISE_MATCH
                for k in range(operation.length):
                    if pairwiseAlignment.contig1 <= pairwiseAlignment.contig2:
                        if pairwiseAlignment.contig1 != pairwiseAlignment.contig2 or i != j: #Avoid self alignments
                            pairsSet.add((pairwiseAlignment.contig1, i, pairwiseAlignment.contig2, j))
                    else:
                        pairsSet.add((pairwiseAlignment.contig2, j, pairwiseAlignment.contig1, i))
                    i += s1
                    j += s2

        if pairwiseAlignment.strand1:
            assert i == pairwiseAlignment.end1
        else:
            assert i == pairwiseAlignment.end1-1

        if pairwiseAlignment.strand2:
            assert j == pairwiseAlignment.end2
        else:
            assert j == pairwiseAlignment.end2-1

        #assert j == pairwiseAlignment.end2
    fileHandle.close()
    return (pairsSet, totalHits)

def runNaiveBlast(seqFile1, seqFile2, outputFile, tempDir, lastzArguments=""):
    """Runs the blast command in a very naive way (not splitting things up).
    """
    startTime = time.time()
    tmpSeqFile1 = os.path.join(tempDir, "seq1.fa")
    tmpSeqFile2 = os.path.join(tempDir, "seq2.fa")
    shutil.copyfile(seqFile1, tmpSeqFile1)
    shutil.copyfile(seqFile2, tmpSeqFile2)
    runLastz(tmpSeqFile1, tmpSeqFile2, alignmentsFile=outputFile, lastzArguments=lastzArguments)
    return time.time()-startTime

def runCactusBlast(sequenceFiles, alignmentsFile, toilDir,
                   chunkSize=None, overlapSize=None,
                   compressFiles=None,
                   lastzMemory=None,
                   targetSequenceFiles=None):

    options = Job.Runner.getDefaultOptions(toilDir)
    options.logLevel = getTestLogLevel()
    blastOptions = BlastOptions(chunkSize=chunkSize, overlapSize=overlapSize,
                                compressFiles=compressFiles,
                                memory=lastzMemory)
    with Toil(options) as toil:
        seqIDs = [toil.importFile(makeURL(seqFile)) for seqFile in sequenceFiles]

        if targetSequenceFiles:
            targetSeqIDs = [toil.importFile(makeURL(seqFile)) for seqFile in targetSequenceFiles]
            rootJob = BlastSequencesAgainstEachOther(sequenceFileIDs1=seqIDs, sequenceFileIDs2=targetSeqIDs, blastOptions=blastOptions)
        else:
            rootJob = BlastSequencesAllAgainstAll(seqIDs, blastOptions)
        alignmentsID = toil.start(rootJob)
        toil.exportFile(alignmentsID, makeURL(alignmentsFile))

def runCactusBlastIngroupsAndOutgroups(ingroups, outgroups, alignmentsFile, toilDir, outgroupFragmentPaths=None, ingroupCoveragePaths=None, chunkSize=250000, overlapSize=10000,
                   compressFiles=None,
                   lastzMemory=None):
    options = Job.Runner.getDefaultOptions(toilDir)
    options.disableCaching = True
    options.logLevel = getTestLogLevel()
    blastOptions = BlastOptions(chunkSize=chunkSize, overlapSize=overlapSize,
                                compressFiles=compressFiles,
                                memory=lastzMemory)
    with Toil(options) as toil:
        ingroupIDs = [toil.importFile(makeURL(ingroup)) for ingroup in ingroups]
        outgroupIDs = [toil.importFile(makeURL(outgroup)) for outgroup in outgroups]
        rootJob = BlastIngroupsAndOutgroups(blastOptions, ingroups, ingroupIDs, outgroups, outgroupIDs)
        blastResults = toil.start(rootJob)
        alignmentsID = blastResults[0]
        toil.exportFile(alignmentsID, makeURL(alignmentsFile))
        outgroupFragmentIDs = blastResults[1]
        ingroupCoverageIDs = blastResults[2]

        if outgroupFragmentPaths:
            assert len(outgroupFragmentIDs) == len(outgroupFragmentPaths)
            for outgroupFragmentID, outgroupFragmentPath in zip(outgroupFragmentIDs, outgroupFragmentPaths):
                toil.exportFile(outgroupFragmentID, makeURL(outgroupFragmentPath))
        if ingroupCoveragePaths:
            assert len(ingroupCoverageIDs) == len(ingroupCoveragePaths)
            for ingroupCoverageID, ingroupCoveragePath in zip(ingroupCoverageIDs, ingroupCoveragePaths):
                toil.exportFile(ingroupCoverageID, makeURL(ingroupCoveragePath))


def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()

if __name__ == '__main__':
    unittest.main()
