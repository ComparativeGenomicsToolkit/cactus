import unittest
import sys
import os
import xml.etree.ElementTree as ET
import random

from sonLib.bioio import logger
from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import system
from sonLib.bioio import getTempDirectory
from sonLib.bioio import fastaRead

from cactus.cactus_common import runCactusSetup
from cactus.cactus_common import getRandomCactusInputs

class TestCase(unittest.TestCase):

    def setUp(self):
        self.testNo = TestStatus.getTestSetup()
        self.tempFiles = []
        self.tempDir = getTempDirectory()
        self.tempReconstructionDirectory = os.path.join(self.tempDir, "tempReconstruction")
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        for tempFile in self.tempFiles:
            os.remove(tempFile)
        unittest.TestCase.tearDown(self)
        system("rm -rf %s" % self.tempDir)
        
    def testCactusSetup(self):
        for test in xrange(self.testNo): 
            ##########################################
            #Make random inputs
            ##########################################
            
            sequenceNumber = random.choice(xrange(100))
            sequenceDirs, newickTreeString = getRandomCactusInputs(tempDir=getTempDirectory(self.tempDir), sequenceNumber=sequenceNumber)
            cactusTempDir=getTempDirectory(self.tempDir)
            
            runCactusSetup(self.tempReconstructionDirectory, sequenceDirs, 
                               newickTreeString, cactusTempDir, debug=True)
            
            runCactusSetup(self.tempReconstructionDirectory, sequenceDirs, 
                               newickTreeString, cactusTempDir, debug=True) #Run it twice to check the job is blockic.
    
            system("rm -rf %s" % self.tempReconstructionDirectory)
            
            logger.info("Finished test of cactus_setup.py")
            

def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
