#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import os
import sys
import random
import time

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
from sonLib.bioio import getLogLevelString
from sonLib.bioio import TestStatus
from sonLib.bioio import catFiles
from sonLib.bioio import popenCatch
from cactus.shared.test import parseCactusSuiteTestOptions
from cactus.shared.common import runCactusBlast
from cactus.blast.cactus_blast import decompressFastaFile, compressFastaFile

from cactus.shared.common import runToilStatusAndFailIfNotComplete

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
        encodeRegions = [ "ENm00" + str(i) for i in xrange(1,2) ] #, 2) ] #Could go to six
        species = ("human", "mouse", "dog")
        #Other species to try "rat", "monodelphis", "macaque", "chimp"
        for encodeRegion in encodeRegions:
            regionPath = os.path.join(self.encodePath, encodeRegion)
            for i in xrange(len(species)):
                species1 = species[i]
                for species2 in species[i+1:]:
                    seqFile1 = os.path.join(regionPath, "%s.%s.fa" % (species1, encodeRegion))
                    seqFile2 = os.path.join(regionPath, "%s.%s.fa" % (species2, encodeRegion))
                    
                    #Run the random
                    runNaiveBlast(seqFile1, seqFile2, self.tempOutputFile)
                    logger.info("Ran the naive blast okay")
                    
                    #Run the blast
                    toilDir = os.path.join(getTempDirectory(self.tempDir), "toil")
                    if blastMode == "allAgainstAll":
                        runCactusBlast([ seqFile1, seqFile2 ], self.tempOutputFile2, toilDir,
                                       chunkSize=500000, overlapSize=10000)
                    else:
                        runCactusBlast([ seqFile1 ], self.tempOutputFile2, toilDir,
                                       chunkSize=500000, overlapSize=10000, targetSequenceFiles=[ seqFile2 ])
                    #runToilStatusAndFailIfNotComplete(toilDir)
                    system("rm -rf %s " % toilDir)    
                    logger.info("Ran cactus_blast okay")
                    logger.critical("Comparing cactus_blast and naive blast; using mode: %s" % blastMode)
                    compareResultsFile(self.tempOutputFile, self.tempOutputFile2)

    def testBlastEncodeAllAgainstAll(self):
        """For each encode region, for set of pairwise species, run 
        cactus_blast.py in all-against-all mode. 
        """
        self.runComparisonOfBlastScriptVsNaiveBlast(blastMode="allAgainstAll")

    def testBlastEncode(self):
        """For each encode region, for set of pairwise species, run 
        cactus_blast.py in one set of sequences against another set mode. 
        """
        self.runComparisonOfBlastScriptVsNaiveBlast(blastMode="againstEachOther")

    def testAddingOutgroupsImprovesResult(self):
        """Run blast on "ingroup" and "outgroup" encode regions, and ensure
        that adding an extra outgroup only adds alignments if
        possible, and doesn't lose any
        """
        encodeRegions = [ "ENm00" + str(i) for i in xrange(1,2) ]
        ingroups = ["human", "macaque"]
        outgroups = ["rabbit", "dog", "rat", "platypus", "xenopus", "fugu"]
        # subselect 4 random ordered outgroups
        outgroups = [outgroups[i] for i in sorted(random.sample(xrange(len(outgroups)), 4))]
        for encodeRegion in encodeRegions:
            regionPath = os.path.join(self.encodePath, encodeRegion)
            ingroupPaths = map(lambda x: os.path.join(regionPath, x + "." + encodeRegion + ".fa"), ingroups)
            outgroupPaths = map(lambda x: os.path.join(regionPath, x + "." + encodeRegion + ".fa"), outgroups)
            results = []
            for numOutgroups in xrange(1,5):
                # Align w/ increasing numbers of outgroups
                subResults = getTempFile()
                subOutgroupPaths = outgroupPaths[:numOutgroups]
                print "aligning %s vs %s" % (",".join(ingroupPaths), ",".join(subOutgroupPaths))
                tmpToil = "./outgroupToil"
                system("cactus_blast.py %s --ingroups %s --outgroups %s --cigars %s" % (tmpToil, ",".join(ingroupPaths), ",".join(subOutgroupPaths), subResults))
                system("rm -fr %s" % (tmpToil))
                results.append(subResults)

            # Print diagnostics about coverage
            for i, subResults in enumerate(results):
                for ingroup, ingroupPath in zip(ingroups, ingroupPaths):
                    coveredBases = popenCatch("cactus_coverage %s %s | awk '{ total += $3 - $2 } END { print total }'" % (ingroupPath, subResults))
                    print "covered bases on %s using %d outgroups: %s" % (ingroup, i + 1, coveredBases)

            resultsSets = map(lambda x : loadResults(x), results)
            for i, moreOutgroupsResults in enumerate(resultsSets[1:]):
                # Make sure the results from (n+1) outgroups are
                # (very nearly) a superset of the results from n outgroups
                print "Using %d addl outgroup(s):" % (i + 1)
                comparator =  ResultComparator(resultsSets[0], moreOutgroupsResults)
                print comparator
                self.assertTrue(comparator.sensitivity >= 0.99)

            # Ensure that the new alignments don't cover more than
            # x% of already existing alignments to human
            for i in xrange(1, len(resultsSets)):
                prevResults = resultsSets[i-1][0]
                curResults = resultsSets[i][0]
                prevResultsHumanPos = set(map(lambda x: (x[0], x[1]) if "human" in x[0] else (x[2], x[3]), filter(lambda x: "human" in x[0] or "human" in x[2], prevResults)))
                newAlignments = curResults.difference(prevResults)
                newAlignmentsHumanPos =  set(map(lambda x: (x[0], x[1]) if "human" in x[0] else (x[2], x[3]), filter(lambda x: "human" in x[0] or "human" in x[2], newAlignments)))
                print "addl outgroup %d:" % i
                print "bases re-covered: %f (%d)" % (len(newAlignmentsHumanPos.intersection(prevResultsHumanPos))/float(len(prevResultsHumanPos)), len(newAlignmentsHumanPos.intersection(prevResultsHumanPos)))
            for subResult in results:
                os.remove(subResult)

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
        outgroupPaths = map(lambda x: os.path.join(regionPath, x + "." + encodeRegion + ".fa"), outgroups)
        # Run in "set against set" mode, aligning the entire ingroup
        # vs each outgroup
        runCactusBlast([ingroupPath], self.tempOutputFile, os.path.join(self.tempDir, "setVsSetToil"),
                       chunkSize=500000, overlapSize=10000,
                       targetSequenceFiles=outgroupPaths)
        # Run in "ingroup vs outgroups" mode, aligning the ingroup vs
        # the outgroups in order, trimming away sequence that's
        # already been aligned.
        system("cactus_blast.py %s/outgroupToil --ingroups %s --outgroups %s --cigars %s" % (self.tempDir, ingroupPath, ",".join(outgroupPaths), self.tempOutputFile2))

        # Get the coverage on the ingroup, in bases, from each run.
        coverageSetVsSet = int(popenCatch("cactus_coverage %s %s | awk '{ total +=  $3 - $2} END { print total }'" % (ingroupPath, self.tempOutputFile)))
        coverageIngroupVsOutgroups = int(popenCatch("cactus_coverage %s %s | awk '{ total +=  $3 - $2} END { print total }'" % (ingroupPath, self.tempOutputFile2)))

        print "total coverage on human (set vs set mode, %d outgroups): %d" % (len(outgroups), coverageSetVsSet)
        print "total coverage on human (ingroup vs outgroup mode, %d outgroups): %d" % (len(outgroups), coverageIngroupVsOutgroups)

        # Make sure we're getting a reasonable fraction of the
        # alignments when using the trimming strategy.
        self.assertTrue(float(coverageIngroupVsOutgroups)/coverageSetVsSet >= 0.95)

        # Get the coverage on the ingroup, in bases, from just the
        # last outgroup. Obviously this should be much higher in set
        # vs set mode than in ingroup vs outgroup mode.
        coverageFromLastOutgroupSetVsSet = int(popenCatch("grep %s %s | cactus_coverage %s /dev/stdin | awk '{ total +=  $3 - $2} END { print total }'" % (outgroups[-1], self.tempOutputFile, ingroupPath)))
        coverageFromLastOutgroupInVsOut = int(popenCatch("grep %s %s | cactus_coverage %s /dev/stdin | awk '{ total +=  $3 - $2} END { print total }'" % (outgroups[-1], self.tempOutputFile2, ingroupPath)))

        print "total coverage on human from last outgroup in set (%s) (set vs set mode): %d" % (outgroups[-1], coverageFromLastOutgroupSetVsSet)
        print "total coverage on human from last outgroup in set (%s) (ingroup vs outgroup mode): %d" % (outgroups[-1], coverageFromLastOutgroupInVsOut)

        self.assertTrue(float(coverageFromLastOutgroupInVsOut)/coverageFromLastOutgroupSetVsSet <= 0.10)
        
    def testBlastParameters(self):
        """Tests if changing parameters of lastz creates results similar to the desired default.
        """
        encodeRegion = "ENm001"
        species = ("human", "mouse", "dog")
        #Other species to try "rat", "monodelphis", "macaque", "chimp"
        regionPath = os.path.join(self.encodePath, encodeRegion)
        for i in xrange(len(species)):
            species1 = species[i]
            for species2 in species[i+1:]:
                seqFile1 = os.path.join(regionPath, "%s.%s.fa" % (species1, encodeRegion))
                seqFile2 = os.path.join(regionPath, "%s.%s.fa" % (species2, encodeRegion))
                
                #Run the random
                runNaiveBlast(seqFile1, seqFile2, self.tempOutputFile2,
                              lastzOptions="--nogapped --hspthresh=3000 --ambiguous=iupac")
                #Run the blast
                runNaiveBlast(seqFile1, seqFile2, self.tempOutputFile, 
                              lastzOptions="--nogapped --step=3 --hspthresh=3000 --ambiguous=iupac")
                
                logger.critical("Comparing blast settings")
                compareResultsFile(self.tempOutputFile, self.tempOutputFile2, 0.7)

    def testBlastRandom(self):
        """Make some sequences, put them in a file, call blast with random parameters 
        and check it runs okay.
        """
        tempSeqFile = os.path.join(self.tempDir, "tempSeq.fa")
        self.tempFiles.append(tempSeqFile)
        for test in xrange(self.testNo):
            seqNo = random.choice(xrange(0, 10))
            seq = getRandomSequence(8000)[1]
            fileHandle = open(tempSeqFile, 'w')
            for fastaHeader, seq in [ (str(i), mutateSequence(seq, 0.3*random.random())) for i in xrange(seqNo) ]:
                if random.random() > 0.5:
                    seq = reverseComplement(seq)
                fastaWrite(fileHandle, fastaHeader, seq)
            fileHandle.close()
            chunkSize = random.choice(xrange(500, 9000))
            overlapSize = random.choice(xrange(2, 100))
            toilDir = os.path.join(getTempDirectory(self.tempDir), "toil")
            runCactusBlast([ tempSeqFile ], self.tempOutputFile, toilDir, chunkSize, overlapSize)
            #runToilStatusAndFailIfNotComplete(toilDir)
            if getLogLevelString() == "DEBUG":
                system("cat %s" % self.tempOutputFile)
            system("rm -rf %s " % toilDir)

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
                for k in xrange(operation.length):
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

def runNaiveBlast(seqFile1, seqFile2, outputFile, 
                  blastString="cactus_lastz --format=cigar OPTIONS SEQ_FILE_1[multiple][nameparse=darkspace] SEQ_FILE_2[nameparse=darkspace] > CIGARS_FILE", 
                  lastzOptions=""):
    """Runs the blast command in a very naive way (not splitting things up).
    """
    open(outputFile, 'w').close() #Ensure is empty of results
    command = blastString.replace("OPTIONS", lastzOptions).replace("CIGARS_FILE", outputFile).replace("SEQ_FILE_1", seqFile1).replace("SEQ_FILE_2", seqFile2)
    startTime = time.time()
    system(command)
    return time.time()-startTime

def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
