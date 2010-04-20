import unittest
import sys
import os
from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import system
from sonLib.bioio import getTempFile
from cactus.shared.cactus_common import getRandomNetDisk

class TestCase(unittest.TestCase):

    def setUp(self):
        self.testNo = TestStatus.getTestSetup()
        self.tempFiles = []
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        for tempFile in self.tempFiles:
            os.remove(tempFile)
        unittest.TestCase.tearDown(self)
        
    def testCactusMafGenerator(self):
        for test in xrange(self.testNo):
            netDisk = getRandomNetDisk(tempDir=self.tempFileDir)
            mafFile = getTempFile()
            runCactusMAFGenerator(netDisk, mafFile)
            

def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()