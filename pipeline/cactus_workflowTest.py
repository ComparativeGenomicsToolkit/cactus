"""Tests the core pipeline.
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
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import getCactusInputs_encode
from cactus.shared.test import getCactusInputs_chromosomeX
from cactus.shared.test import runWorkflow_TestScript

class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.testNo = TestStatus.getTestSetup(10, 1, 1, 1)
        self.sequenceNumber = TestStatus.getTestSetup(5, 50, 50, 50)
        if os.system("parasol status") == 0:
            self.batchSystem = "parasol"
        else:
            self.batchSystem = "single_machine"
        unittest.TestCase.setUp(self)
        
    def testCactusWorkflow_Random(self):
        """Runs the tests across some simulated regions.
        """
        tempDir = getTempDirectory(os.getcwd())
        for test in xrange(self.testNo): 
            sequences, newickTreeString = getCactusInputs_random(tempDir)
            runWorkflow_TestScript(sequences, newickTreeString, tempDir, 
                                   batchSystem=self.batchSystem,)
            logger.info("Finished random test %i" % test)
        system("rm -rf %s" % tempDir)

    def testCactusWorkflow_Blanchette(self): 
        """Runs the workflow on blanchette's simulated (colinear) regions.
        """
        if TestStatus.getTestStatus() in (TestStatus.TEST_MEDIUM,):
            tempDir = getTempDirectory(os.getcwd())
            for region in xrange(0, 5):
                sequences, newickTreeString = getCactusInputs_blanchette(region)
                outputDir = os.path.join(TestStatus.getPathToDataSets(), "cactus", 
                                 "blanchettesRegionsTest", str(region))
                runWorkflow_TestScript(sequences, newickTreeString, tempDir,
                                       outputDir=outputDir,
                                       batchSystem=self.batchSystem,
                                       buildCactusPDF=True, buildAdjacencyPDF=True,
                                       makeCactusTreeStats=True)
                logger.info("Finished blanchette test %i" % region)
            system("rm -rf %s" % tempDir)
                
    def testCactusWorkflow_Encode(self): 
        """Run the workflow on the encode pilot regions.
        """
        if TestStatus.getTestStatus() in (TestStatus.TEST_LONG,):
            tempDir = getTempDirectory(os.getcwd())
            for encodeRegion in [ "ENm00" + str(i) for i in xrange(1, 5) ]:
                sequences, newickTreeString = getCactusInputs_encode(encodeRegion)
                outputDir = os.path.join(TestStatus.getPathToDataSets(), "cactus", "encodeRegionsTest", encodeRegion)
                runWorkflow_TestScript(sequences, newickTreeString, tempDir,
                                       batchSystem=self.batchSystem,
                                       outputDir=outputDir, makeCactusTreeStats=True)
                logger.info("Finished encode test %i" % encodeRegion)
            system("rm -rf %s" % tempDir)
    
    def testCactusWorkflow_Chromosomes(self):
        """Run the workflow on mammalian chromsome X
        """
        if TestStatus.getTestStatus() in (TestStatus.TEST_VERY_LONG,):
            tempDir = getTempDirectory(os.getcwd())
            sequences, newickTreeString = getCactusInputs_chromosomeX()
            outputDir = os.path.join(TestStatus.getPathToDataSets(), "cactus", "chrX")
            runWorkflow_TestScript(sequences, newickTreeString, tempDir,
                                   batchSystem=self.batchSystem,
                                   outputDir=outputDir, makeCactusTreeStats=True)
            system("rm -rf %s" % tempDir)
            logger.info("Finished chromsome test")
    
def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
