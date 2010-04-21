import unittest
import sys
import os
import random

from sonLib.bioio import logger
from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import system
from sonLib.bioio import getTempDirectory

from cactus.shared.common import runCactusSetup
from cactus.shared.test import getCactusInputs_random

class TestCase(unittest.TestCase):

    def setUp(self):
        self.testNo = TestStatus.getTestSetup()
        unittest.TestCase.setUp(self)
       
    def testCactusSetup(self):
        """Creates a bunch of random inputs and then passes them to cactus setup.
        """
        for test in xrange(self.testNo): 
            tempDir = getTempDirectory(os.getcwd())
            sequenceNumber = random.choice(xrange(100))
            sequences, newickTreeString = getCactusInputs_random(tempDir=tempDir, sequenceNumber=sequenceNumber)
            netDisk = os.path.join(tempDir, "netDisk")
            runCactusSetup(netDisk, sequences, newickTreeString, debug=True)
            runCactusSetup(netDisk, sequences, newickTreeString, debug=True)
            system("rm -rf %s" % tempDir)
            logger.info("Finished test %i of cactus_setup.py", test)

def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
