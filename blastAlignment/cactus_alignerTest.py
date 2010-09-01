import unittest
import sys
import os
from sonLib.bioio import logger
from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import system
from sonLib.bioio import getTempDirectory

from cactus.shared.test import getCactusInputs_random
from cactus.shared.config import CactusWorkflowExperiment

from cactus.shared.common import runCactusSetup
from cactus.shared.common import runCactusAligner

class TestCase(unittest.TestCase):

    def setUp(self):
        self.testNo = TestStatus.getTestSetup(1,5, 50, 0)
        unittest.TestCase.setUp(self)
       
    def testCactusAligner(self):
        """Tests cactus aligner, using the dummy aligner.
        """
        for test in xrange(self.testNo): 
            runAligner()
            
    def testCactusAlignerWithBatchBlast(self):
        """Tests cactus aligner, using the dummy aligner.
        """
        for test in xrange(self.testNo): 
            runAligner(useDummy=False)
    
def runAligner(useDummy=True):       
    ##########################################
    #Make inputs
    ##########################################
    
    tempDir = getTempDirectory(".")
    tempAlignmentFile = os.path.join(tempDir, "tempAlignment.cigar")
    
    sequenceDirs, newickTreeString = getCactusInputs_random(tempDir=getTempDirectory(tempDir))
    outputDir = getTempDirectory(tempDir)
    
    experiment = CactusWorkflowExperiment(sequenceDirs, newickTreeString, outputDir)
    cactusDiskDatabaseString = experiment.getDatabaseString()
    
    runCactusSetup(cactusDiskDatabaseString, sequenceDirs, newickTreeString)
    logger.info("Setup the test, running cactus aligner")
    
    ##########################################
    #Do actual run
    ##########################################
    
    runCactusAligner(cactusDiskDatabaseString, tempAlignmentFile, 
                     tempDir=getTempDirectory(tempDir), useDummy=useDummy)
    runCactusAligner(cactusDiskDatabaseString, tempAlignmentFile, \
                     tempDir=getTempDirectory(tempDir), useDummy=useDummy) #Run it twice to check blockicity
    logger.info("Ran cactus aligner")
    
    system("cat %s" % tempAlignmentFile)
    
    #Cleanup
    experiment.cleanupDatabase()
    system("rm -rf %s" % tempDir)

def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
