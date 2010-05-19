import unittest
import sys
import os
from sonLib.bioio import logger
from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import system
from sonLib.bioio import getTempDirectory
from sonLib.bioio import getTempFile

from cactus.shared.test import getCactusInputs_random

from cactus.shared.common import runCactusSetup
from cactus.shared.common import runCactusAligner

from workflow.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete 

class TestCase(unittest.TestCase):

    def setUp(self):
        self.testNo = TestStatus.getTestSetup(1,5, 50, 0)
        self.tempFiles = []
        self.tempDir = getTempDirectory(".")
        self.tempReconstructionDirectory = os.path.join(self.tempDir, "tempReconstruction")
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        for tempFile in self.tempFiles:
            os.remove(tempFile)
        unittest.TestCase.tearDown(self)
        system("rm -rf %s" % self.tempDir)
        
    def testCactusAligner(self):
        """Tests cactus aligner, using the dummy aligner.
        """
        tempAlignmentFile = getTempFile()
        self.tempFiles.append(tempAlignmentFile)
        for test in xrange(self.testNo): 
            runAligner(self.tempDir, self.tempReconstructionDirectory, tempAlignmentFile)
            
    def testCactusAlignerWithBatchBlast(self):
        """Tests cactus aligner, using the dummy aligner.
        """
        tempAlignmentFile = getTempFile()
        self.tempFiles.append(tempAlignmentFile)
        for test in xrange(self.testNo): 
            runAligner(self.tempDir, self.tempReconstructionDirectory, tempAlignmentFile, useDummy=False)
    
def runAligner(tempDir, tempReconstructionDir, tempAlignmentFile, useDummy=True):       
    ##########################################
    #Make inputs
    ##########################################
    
    sequenceDirs, newickTreeString = getCactusInputs_random(tempDir=getTempDirectory(tempDir))
    runCactusSetup(tempReconstructionDir, sequenceDirs, newickTreeString)
    logger.info("Setup the test, running cactus aligner")
    
    ##########################################
    #Do actual run
    ##########################################
    
    runCactusAligner(tempReconstructionDir, tempAlignmentFile, 
                     tempDir=getTempDirectory(tempDir), useDummy=useDummy)
    runCactusAligner(tempReconstructionDir, tempAlignmentFile, \
                     tempDir=getTempDirectory(tempDir), useDummy=useDummy) #Run it twice to check blockicity
    logger.info("Ran cactus aligner")
    
    system("cat %s" % tempAlignmentFile)

def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
