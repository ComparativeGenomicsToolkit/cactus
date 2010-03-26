import unittest
import os
import sys

from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger
from sonLib.bioio import system

from cactus.shared.cactus_common import getRandomCactusInputs
from cactus.shared.cactus_common import runCactusWorkflow
from cactus.shared.cactus_common import runCactusCheck

from workflow.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete

class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.testNo = TestStatus.getTestSetup(1, 5, 100, 100)
        self.sequenceNumber = TestStatus.getTestSetup(5, 50, 50, 100)
        self.tempFiles = []
        self.tempDir = getTempDirectory(os.getcwd())
        self.tempReconstructionDirectory = os.path.join(self.tempDir, "tempReconstruction")
        if os.system("parasol status") == 0:
            self.batchSystem = "parasol"
        else:
            self.batchSystem = "single_machine"
        #self.batchSystem = "single_machine"
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        for tempFile in self.tempFiles:
            os.remove(tempFile)
        unittest.TestCase.tearDown(self)
        system("rm -rf %s" % self.tempDir)
        
    def testCactusReference_Random(self):
        """Creates a reference genome for random regions.
        """
        for test in xrange(self.testNo): 
            sequenceDirs, newickTreeString = getRandomCactusInputs(tempDir=getTempDirectory(self.tempDir), sequenceNumber=self.sequenceNumber)
            jobTreeDir = os.path.join(getTempDirectory(self.tempDir), "jobTree")
            runCactusWorkflow(self.tempReconstructionDirectory, sequenceDirs, newickTreeString, jobTreeDir, logLevel='INFO',
                               batchSystem=self.batchSystem, buildReference=True)
            #runCactusCheck(self.tempReconstructionDirectory, checkReference=True)
            
            #Run the checker to check the file is okay.
            runJobTreeStatusAndFailIfNotComplete(jobTreeDir)
            #Cleanup this test
            system("rm -rf %s %s" % (self.tempReconstructionDirectory, jobTreeDir))
            logger.info("Finished this round of test")

def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
