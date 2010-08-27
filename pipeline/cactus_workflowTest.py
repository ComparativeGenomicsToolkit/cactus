"""Tests the core pipeline.
"""

import unittest
import os
import sys

from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import getCactusInputs_encode
from cactus.shared.test import getCactusInputs_chromosomeX
from cactus.shared.test import runWorkflow_multipleExamples

class TestCase(unittest.TestCase):
    
    def setUp(self):
        if os.system("parasol status") == 0:
            self.batchSystem = "parasol"
        else:
            self.batchSystem = "single_machine"
        unittest.TestCase.setUp(self)
        
    def testCactus_Random(self):
        runWorkflow_multipleExamples(getCactusInputs_random, 
                                     testNumber=TestStatus.getTestSetup(),
                                     batchSystem=self.batchSystem)
        
    def testCactus_Blanchette(self):
        outputDir = os.path.join(TestStatus.getPathToDataSets(), "cactus", 
                                 "blanchettesRegionsTest")
        runWorkflow_multipleExamples(getCactusInputs_blanchette, 
                                     outputDir=outputDir,
                                     testNumber=5,
                                     testRestrictions=(TestStatus.TEST_MEDIUM,),
                                     batchSystem=self.batchSystem,
                                     buildCactusPDF=True,
                                     makeCactusTreeStats=True, makeMAFs=True, buildReferencePDF=True)
                
    def testCactus_Encode(self): 
        outputDir = os.path.join(TestStatus.getPathToDataSets(), "cactus", "encodeRegionsTest")
        runWorkflow_multipleExamples(getCactusInputs_encode, 
                                     outputDir=outputDir,
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_LONG,),
                                     batchSystem=self.batchSystem,
                                     makeCactusTreeStats=True, makeMAFs=True)
    
    def testCactus_Chromosomes(self):
        outputDir = os.path.join(TestStatus.getPathToDataSets(), "cactus", "chrX")
        runWorkflow_multipleExamples(getCactusInputs_chromosomeX, 
                                     outputDir=outputDir,
                                     testRestrictions=(TestStatus.TEST_VERY_LONG,),
                                     batchSystem=self.batchSystem,
                                     makeCactusTreeStats=True, makeMAFs=True)
    
def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
