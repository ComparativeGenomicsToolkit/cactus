#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
"""Tests the core pipeline.
"""

import unittest
import os
import sys

from cactus.shared.test import parseCactusSuiteTestOptions
from sonLib.bioio import TestStatus

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import getCactusInputs_encode
from cactus.shared.test import getCactusInputs_chromosomeX

from cactus.shared.test import runWorkflow_multipleExamples

from cactus.shared.test import getBatchSystem

class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.batchSystem = "singleMachine"
        if getBatchSystem() != None:
            self.batchSystem = getBatchSystem()
        unittest.TestCase.setUp(self)
        
    def testCactus_Random(self):
        runWorkflow_multipleExamples(getCactusInputs_random, 
                                     testNumber=2,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True)
        
    def testCactus_Blanchette(self):
        outputDir = os.path.join(TestStatus.getPathToDataSets(), "cactus", 
                                 "blanchettesRegionsTest")
        runWorkflow_multipleExamples(getCactusInputs_blanchette, 
                                     outputDir=outputDir,
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_MEDIUM,),
                                     batchSystem=self.batchSystem,
                                     buildCactusPDF=True,
                                     makeCactusTreeStats=True, makeMAFs=True, buildReferencePDF=True, buildJobTreeStats=True)
                
    def testCactus_Encode(self): 
        outputDir = os.path.join(TestStatus.getPathToDataSets(), "cactus", "encodeRegionsTest")
        runWorkflow_multipleExamples(getCactusInputs_encode, 
                                     outputDir=outputDir,
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_LONG,),
                                     batchSystem=self.batchSystem,
                                     makeCactusTreeStats=True, makeMAFs=True, buildJobTreeStats=True)
    
    def testCactus_Chromosomes(self):
        outputDir = os.path.join(TestStatus.getPathToDataSets(), "cactus", "chrX")
        runWorkflow_multipleExamples(getCactusInputs_chromosomeX, 
                                     outputDir=outputDir,
                                     testRestrictions=(TestStatus.TEST_VERY_LONG,),
                                     batchSystem=self.batchSystem,
                                     makeCactusTreeStats=True, makeMAFs=True)
    
def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
