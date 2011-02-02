import unittest
import os
import sys
import random

from sonLib.bioio import system

from sonLib.bioio import logger
from sonLib.bioio import fastaWrite
from sonLib.bioio import getRandomSequence
from sonLib.bioio import mutateSequence
from sonLib.bioio import reverseComplement

from sonLib.bioio import getTempFile
from sonLib.bioio import getTempDirectory

from sonLib.bioio import TestStatus
from cactus.shared.test import parseCactusSuiteTestOptions

from sonLib.bioio import cigarRead
from sonLib.bioio import PairwiseAlignment

from workflow.jobTree.test.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete

class TestCase(unittest.TestCase):
    def setUp(self):
        self.testNo = TestStatus.getTestSetup(1, 5, 10, 100)
        self.tempDir = getTempDirectory(os.getcwd())
        self.tempFiles = []
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        for tempFile in self.tempFiles:
            os.remove(tempFile)
        unittest.TestCase.tearDown(self)
        system("rm -rf %s" % self.tempDir)
    
    def testChromosomes(self):
        """Tests cactus batch on the alignment 4 whole chromosome X's, human, chimp, mouse, dog.
        """
        if TestStatus.getTestStatus() in (TestStatus.TEST_VERY_LONG,):
            tempOutputFile = getTempFile()
            self.tempFiles.append(tempOutputFile)
            jobTreeDir = os.path.join(getTempDirectory(self.tempDir), "jobTree")
            chrXPath = os.path.join(TestStatus.getPathToDataSets(), "chr_x")
            seqs = ("hg18.fa", "panTro2.fa", "dog_chrX.fa", "mouse_chrX.fa")
            options = makeTopLevelBlastOptions()
            runCactusBatch([ os.path.join(chrXPath, seqFile) for seqFile in seqs ], tempOutputFile, jobTreeDir,
                   chunkSize=options.chunkSize, overlapSize=options.overlapSize, chunksPerJob=options.chunksPerJob,
                   logLevel="DEBUG", retryCount=0, batchSystem="single_machine", 
                   blastString=options.blastString, selfBlastString=options.selfBlastString)
        
    def testBatchEncode(self):
        """For each encode region, for set of pairwise species, run 
        cactus_batch.py. We compare the output with a naive run of the blast program, to check the results are nearly
        equivalent.
        """
        tempOutputFile = getTempFile()
        self.tempFiles.append(tempOutputFile)
        tempOutputFile2 = getTempFile()
        self.tempFiles.append(tempOutputFile2) 
        
        if TestStatus.getTestStatus() in (TestStatus.TEST_LONG, TestStatus.TEST_VERY_LONG):
            encodePath = os.path.join(TestStatus.getPathToDataSets(), "MAY-2005")
            encodeRegions = [ "ENm00" + str(i) for i in xrange(1, 2) ] #Could go to six
            species = ("human", "mouse", "dog")
            #Other species to try "rat", "monodelphis", "macaque", "chimp"
            for encodeRegion in encodeRegions:
                regionPath = os.path.join(encodePath, encodeRegion)
                for i in xrange(len(species)):
                    species1 = species[i]
                    for species2 in species[i+1:]:
                        seqFile1 = os.path.join(regionPath, "%s.%s.fa" % (species1, encodeRegion))
                        seqFile2 = os.path.join(regionPath, "%s.%s.fa" % (species2, encodeRegion))
                        
                        #Run the random
                        runNaiveBlast([ seqFile1, seqFile2 ], tempOutputFile2, self.tempDir)
                        logger.info("Ran the naive blast okay")
                        
                        results2 = loadResults(tempOutputFile2)
                        logger.info("Loaded naive blast results")
                        
                        checkResultsAreApproximatelyEqual(ResultComparator(results2, results2)) #Dummy check
                        
                        #Run the batch
                        jobTreeDir = os.path.join(getTempDirectory(self.tempDir), "jobTree")
                        runCactusBatch([ seqFile1, seqFile2 ], tempOutputFile, jobTreeDir,
                                       chunkSize=500000, overlapSize=10000)
                        runJobTreeStatusAndFailIfNotComplete(jobTreeDir)
                        system("rm -rf %s " % jobTreeDir)    
                        logger.info("Ran cactus_batch okay")
                        
                        #Now compare the results
                        results1 = loadResults(tempOutputFile)
                        logger.info("Loaded cactus_batch results")
                        
                        checkResultsAreApproximatelyEqual(ResultComparator(results2, results1))
                        logger.info("Compared the naive and batch results, using the naive results as the 'true' results, and the batch results as the predicted results")
    
    def testBatchRandom(self):
        """Make some sequences, put them in a file, call batch with random parameters 
        and check it runs okay.
        """
        tempSeqFile = getTempFile()
        self.tempFiles.append(tempSeqFile)
            
        tempOutputFile = getTempFile()
        self.tempFiles.append(tempOutputFile)
        
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
            chunksPerJob = random.choice(xrange(1, 5))
            jobTreeDir = os.path.join(getTempDirectory(self.tempDir), "jobTree")
            runCactusBatch([ tempSeqFile ], tempOutputFile, jobTreeDir, chunkSize, overlapSize, chunksPerJob)
            runJobTreeStatusAndFailIfNotComplete(jobTreeDir)
            system("cat %s" % tempOutputFile)
            system("rm -rf %s " % jobTreeDir)

def checkResultsAreApproximatelyEqual(resultsComparator):
    """Checks that the comparisons show the two blasts are suitably close together.
    """
    logger.info("Results are: %s" % resultsComparator)
    assert resultsComparator.sensitivity >= 0.95
    assert resultsComparator.specificity >= 0.95
    
class ResultComparator:
    def __init__(self, trueResults, predictedResults):
        """Compares two sets of results and returns a set of statistics comparing them.
        """
        #Totals
        self.trueLength = len(trueResults)
        self.predictedLength = len(predictedResults)
        self.unionSize = len(trueResults.union(predictedResults))
        self.intersectionSize = len(trueResults.intersection(predictedResults))
        self.symmDiff = self.intersectionSize / self.unionSize
        #Sensitivity
        self.trueDifference = float(len(trueResults.difference(predictedResults)))
        self.sensitivity = 1.0 - float(self.trueDifference) / len(trueResults)
        #Specificity
        self.predictedDifference = float(len(predictedResults.difference(trueResults)))
        self.specificity = 1.0 - float(self.predictedDifference) / len(predictedResults)

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
    for pairwiseAlignment in cigarRead(fileHandle):
        i = pairwiseAlignment.start1
        if not pairwiseAlignment.strand1:
            i -= 1
            
        j = pairwiseAlignment.start2
        if not pairwiseAlignment.strand2:
            j -= 1
        
        for operation in pairwiseAlignment.operationList:
            if operation.type == PairwiseAlignment.PAIRWISE_INDEL_X:
                if pairwiseAlignment.strand1:
                    i += operation.length
                else:
                    i -= operation.length
                    
            elif operation.type == PairwiseAlignment.PAIRWISE_INDEL_Y:
                if pairwiseAlignment.strand2:
                    j += operation.length
                else:
                    j -= operation.length
            else:
                assert operation.type == PairwiseAlignment.PAIRWISE_MATCH
                for k in xrange(operation.length):
                    pairsSet.add((pairwiseAlignment.contig1, i, pairwiseAlignment.contig2, j))
                    pairsSet.add((pairwiseAlignment.contig2, j, pairwiseAlignment.contig1, i)) #Add them symmetrically
                    if pairwiseAlignment.strand1:
                        i += 1
                    else:
                        i -= 1
                    if pairwiseAlignment.strand2:
                        j += 1
                    else:
                        j -= 1
        
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
    return pairsSet

def runNaiveBlast(sequenceFiles, outputFile, tempDir,
                  blastString="lastz --format=cigar SEQ_FILE_1[multiple][nameparse=darkspace] SEQ_FILE_2[nameparse=darkspace] > CIGARS_FILE", 
                  selfBlastString="lastz --format=cigar SEQ_FILE[nameparse=darkspace] --self > CIGARS_FILE"):
    """Runs the blast command in a very naive way (not splitting things up).
    """
    open(outputFile, 'w').close() #Ensure is empty of results
    tempResultsFile = getTempFile(suffix=".results", rootDir=tempDir)
    for i in xrange(len(sequenceFiles)):
        seqFile1 = sequenceFiles[i]
        command = selfBlastString.replace("CIGARS_FILE", tempResultsFile).replace("SEQ_FILE", seqFile1)
        system(command)
        logger.info("Ran the self blast okay for sequence: %s" % seqFile1)
        system("cat %s >> %s" % (tempResultsFile, outputFile))
        for j in xrange(i+1, len(sequenceFiles)):
            seqFile2 = sequenceFiles[j]
            command = blastString.replace("CIGARS_FILE", tempResultsFile).replace("SEQ_FILE_1", seqFile1).replace("SEQ_FILE_2", seqFile2)
            system(command)
            logger.info("Ran the blast okay for sequences: %s %s" % (seqFile1, seqFile2))
            system("cat %s >> %s" % (tempResultsFile, outputFile))
    os.remove(tempResultsFile)

def runCactusBatch(sequenceFiles, outputFile, jobTreeDir,
                   chunkSize=1000000, overlapSize=1000, chunksPerJob=1,
                   logLevel="DEBUG", retryCount=0, batchSystem="single_machine",
                   blastString="lastz --format=cigar SEQ_FILE_1[multiple][nameparse=darkspace] SEQ_FILE_2[nameparse=darkspace] > CIGARS_FILE", 
                   selfBlastString="lastz --format=cigar SEQ_FILE[nameparse=darkspace] --self > CIGARS_FILE"):
    command = "cactus_batch.py %s  \
--cigars %s --chunkSize %s --overlapSize %s --chunksPerJob %s --blastString '%s' --selfBlastString '%s' --jobTree %s --logLevel %s \
--retryCount %i --batchSystem %s" % \
            (" ".join(sequenceFiles), outputFile,
             chunkSize, overlapSize, chunksPerJob, blastString, selfBlastString, jobTreeDir, logLevel, retryCount, batchSystem)
    logger.info("Running command : %s" % command)
    system(command)
    logger.info("Ran the cactus_batch command okay")

def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
