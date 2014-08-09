import unittest
import os, sys
from cactus.shared.test import parseCactusSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import logger
from sonLib.bioio import getTempFile, getTempDirectory, system
from cactus.bar.cactus_realignTest import seqFilePairGenerator
from cactus.shared.common import runCactusExpectationMaximisation
from cactus.bar.cactus_expectationMaximisation import Hmm
from sonLib.bioio import getLogLevelString

class TestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)

    def testCactusExpectationMaximisation(self):
        """Runs cactus expectation-maximisation. 
        """
        trial = 0
        for modelType in ("fiveState", "threeState", "threeStateAsymmetric"):
            for seqFile1, seqFile2 in seqFilePairGenerator():
                tempDir = getTempDirectory(rootDir=os.getcwd())
                jobTreeDir = os.path.join(tempDir, "jobTree")
                alignmentsFile = os.path.join(tempDir, "alignments.cigars")
                computeAlignments(seqFile1, seqFile2, alignmentsFile)
                logger.info("Computed alignments for seqs %s and %s" % (seqFile1, seqFile2))
                outputModelFile = os.path.join(tempDir, "outputModel.txt")
                #First run the script to generate a model and do one iteration of EM to 
                #get the likelihood to compare with the final likelihood
                runCactusExpectationMaximisation(sequenceFiles=[ seqFile1, seqFile2 ], 
                                                 alignmentsFile=alignmentsFile, outputModelFile=outputModelFile, 
                                                 modelType=modelType,
                                                 jobTreeDir=jobTreeDir,
                                                 iterations=1, trials=1, randomStart=False, logLevel=getLogLevelString(),
                                                 optionsToRealign="--diagonalExpansion=5 --splitMatrixBiggerThanThis=3000")
                hmm = Hmm.loadHmm(outputModelFile)
                system("rm -rf %s" % jobTreeDir) #Cleanup the old jobTree
                logger.info("For trial %s the likelihood after 1 iteration of EM is %s" % (trial, hmm.likelihood))
                iterations = 5
                runCactusExpectationMaximisation(sequenceFiles=[ seqFile1, seqFile2 ], 
                                                 alignmentsFile=alignmentsFile, outputModelFile=outputModelFile, jobTreeDir=jobTreeDir,
                                                 optionsToRealign="--diagonalExpansion=5 --splitMatrixBiggerThanThis=100",
                                                 iterations=iterations, inputModelFile=outputModelFile, logLevel=getLogLevelString())
                hmm2 = Hmm.loadHmm(outputModelFile)
                logger.info("For trial %s the likelihood after a further %s iterations of EM is %s" % (trial, iterations, hmm2.likelihood))
                self.assertTrue(hmm.likelihood < hmm2.likelihood)
                hmm2.normalise()
                logger.info("Final transitions: %s" % " ".join(map(str, hmm2.transitions)))
                system("rm -rf %s" % tempDir)
                trial += 1
    
    def testCactusExpectationMaximisationMultipleTrials(self):
        """Runs cactus expectation-maximisation with multiple different trials.
        """
        for seqFile1, seqFile2 in seqFilePairGenerator():
            tempDir = getTempDirectory(rootDir=os.getcwd())
            jobTreeDir = os.path.join(tempDir, "jobTree")
            alignmentsFile = os.path.join(tempDir, "alignments.cigars")
            computeAlignments(seqFile1, seqFile2, alignmentsFile)
            logger.info("Computed alignments for seqs %s and %s" % (seqFile1, seqFile2))
            outputModelFile = os.path.join(tempDir, "outputModel.txt")
            #First run the script to generate a model and do one iteration of EM to 
            #get the likelihood to compare with the final likelihood
            runCactusExpectationMaximisation(sequenceFiles=[ seqFile1, seqFile2 ], 
                                             alignmentsFile=alignmentsFile, outputModelFile=outputModelFile, 
                                             jobTreeDir=jobTreeDir,
                                             trials=3,
                                             iterations=5, randomStart=True, logLevel=getLogLevelString(),
                                             optionsToRealign="--diagonalExpansion=5 --splitMatrixBiggerThanThis=100")
            hmm = Hmm.loadHmm(outputModelFile)
            logger.info("After multiple trials and iterations of EM the best likelihood found was %s" % hmm.likelihood)
            system("rm -rf %s" % tempDir)

def computeAlignments(seqFile1, seqFile2, alignmentsFile, lastzArguments="--ambiguous=iupac"):  
    system("cactus_lastz --format=cigar %s %s[multiple][nameparse=darkspace] %s[nameparse=darkspace] > %s" % (lastzArguments, seqFile1, seqFile2, alignmentsFile))

def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
