#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
"""Tests the core pipeline with Evolver simulated phylogenies.
"""

import unittest
import os
import sys
from cactus.shared.test import parseCactusSuiteTestOptions
from sonLib.bioio import TestStatus
from cactus.shared.test import getInputs
from cactus.shared.test import runWorkflow_multipleExamples
from cactus.shared.test import getBatchSystem

class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.batchSystem = "singleMachine"
        if getBatchSystem() != None:
            self.batchSystem = getBatchSystem()
        unittest.TestCase.setUp(self)
        
    def testEvolver_Primates_Small(self):
        inputDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "primates", "small")
        outputDir = os.path.join(TestStatus.getPathToDataSets(), "cactus", "evolver", "primates", "small")
        primateSequences = ("simChimp.fa", "simGorilla.fa", "simHuman.fa", "simOrang.fa")
        runWorkflow_multipleExamples(lambda regionNumber=0, tempDir=None : getInputs(inputDir, primateSequences),
                                     outputDir=outputDir,
                                     testRestrictions=(TestStatus.TEST_MEDIUM,),
                                     batchSystem=self.batchSystem,
                                     buildTrees=False, buildReference=False,
                                     buildJobTreeStats=True)
    
    def testEvolver_Primates_Loci1(self):
        inputDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "primates", "loci1")
        outputDir = os.path.join(TestStatus.getPathToDataSets(), "cactus", "evolver", "primates", "loci1")
        primateSequences = ("simChimp.chr6", "simGorilla.chr6", "simHuman.chr6", "simOrang.chr6")
        runWorkflow_multipleExamples(lambda regionNumber=0, tempDir=None : getInputs(inputDir, primateSequences),
                                     outputDir=outputDir,
                                     testRestrictions=(TestStatus.TEST_LONG,),
                                     batchSystem=self.batchSystem,
                                     buildTrees=False, buildReference=False,
                                     buildJobTreeStats=True)
    
    def testEvolver_Mammals_Loci1(self):
        inputDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "mammals", "loci1")
        outputDir = os.path.join(TestStatus.getPathToDataSets(), "cactus", "evolver", "mammals", "loci1")
        mammalSequences = ("simCow.chr6", "simDog.chr6", "simHuman.chr6", "simMouse.chr6", "simRat.chr6")
        runWorkflow_multipleExamples(lambda regionNumber=0, tempDir=None : getInputs(inputDir, mammalSequences),
                                     outputDir=outputDir,
                                     testRestrictions=(TestStatus.TEST_VERY_LONG,),
                                     batchSystem=self.batchSystem,
                                     buildTrees=False, buildReference=False,
                                     buildJobTreeStats=True)

def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
