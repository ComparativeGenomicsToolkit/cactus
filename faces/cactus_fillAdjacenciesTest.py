"""Tests the adjacency building of the cactus pipeline.
"""

import unittest
import os
import sys

from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger
from sonLib.bioio import system

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import runWorkflow_TestScript

class TestCase(unittest.TestCase):
    def setUp(self):
        self.testNo = TestStatus.getTestSetup(1, 5, 20, 100)
        unittest.TestCase.setUp(self)
    
    def testCactusFillAdjacencies_Random(self):
        """Runs the tests across some simulated regions.
        """
        tempDir = getTempDirectory(os.getcwd())
        for test in xrange(self.testNo): 
            sequences, newickTreeString = getCactusInputs_random(tempDir)
            runWorkflow_TestScript(sequences, newickTreeString, tempDir, 
                                   buildTrees=True, buildReference=False, 
                                   buildFaces=True)
            logger.info("Finished random test %i" % test)
        system("rm -rf %s" % tempDir)
        
def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
