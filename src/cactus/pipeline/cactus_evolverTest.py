#!/usr/bin/env python3

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
"""Tests the core pipeline with Evolver simulated phylogenies.
"""

import unittest
import os
import sys
from sonLib.bioio import TestStatus
from cactus.shared.test import getInputs
from cactus.shared.test import runWorkflow_multipleExamples
from cactus.shared.test import getBatchSystem

@TestStatus.needsTestData
class TestCase(unittest.TestCase):
    def setUp(self):
        self.batchSystem = "singleMachine"
        if getBatchSystem() != None:
            self.batchSystem = getBatchSystem()
        unittest.TestCase.setUp(self)

    @unittest.skip("test was never updated when changes were made to the way ancestors work (ERROR: Couldn't find reference event reference)")
    @TestStatus.needsTestData
    @TestStatus.longLength
    def testEvolver_Primates_Loci1(self):
        inputDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "primates", "loci1")
        primateSequences = ("simChimp.chr6", "simGorilla.chr6", "simHuman.chr6", "simOrang.chr6")
        runWorkflow_multipleExamples(self.id(),
                                     lambda regionNumber=0, tempDir=None : getInputs(inputDir, primateSequences),
                                     batchSystem=self.batchSystem,
                                     buildToilStats=True)

    @unittest.skip("test was never updated when changes were made to the way ancestors work (ERROR: Couldn't find reference event reference)")
    @TestStatus.needsTestData
    @TestStatus.longLength
    def testEvolver_Mammals_Loci1(self):
        inputDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "mammals", "loci1")
        mammalSequences = ("simCow.chr6", "simDog.chr6", "simHuman.chr6", "simMouse.chr6", "simRat.chr6")
        runWorkflow_multipleExamples(self.id(),
                                     lambda regionNumber=0, tempDir=None : getInputs(inputDir, mammalSequences),
                                     batchSystem=self.batchSystem,
                                     buildToilStats=True)

    @unittest.skip("needs missing cactusTestData/evolver/mammals/medium")
    @TestStatus.needsTestData
    @TestStatus.longLength
    def testEvolver_Mammals_Medium(self):
        inputDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "mammals", "medium")
        mammalSequences = ("simCow.masked.fa", "simDog.masked.fa", "simHuman.masked.fa", "simMouse.masked.fa", "simRat.masked.fa")
        runWorkflow_multipleExamples(self.id(),
                                     lambda regionNumber=0, tempDir=None : getInputs(inputDir, mammalSequences),
                                     batchSystem=self.batchSystem,
                                     buildToilStats=True)

    @unittest.skip("needs missing cactusTestData/evolver/primates/large")
    @TestStatus.needsTestData
    @TestStatus.veryLongLength
    def testEvolver_Primates_Large(self):
        inputDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "primates", "large")
        primateSequences = ("simChimp.masked.fa", "simGorilla.masked.fa", "simHuman.masked.fa", "simOrang.masked.fa")
        runWorkflow_multipleExamples(self.id(),
                                     lambda regionNumber=0, tempDir=None : getInputs(inputDir, primateSequences),
                                     batchSystem=self.batchSystem,
                                     buildToilStats=True)

    @unittest.skip("needs missing cactusTestData/evolver/mammals/large")
    @TestStatus.needsTestData
    @TestStatus.veryLongLength
    def testEvolver_Mammals_Large(self):
        inputDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "mammals", "large")
        mammalSequences = ("simCow.masked.fa", "simDog.masked.fa", "simHuman.masked.fa", "simMouse.masked.fa", "simRat.masked.fa")
        runWorkflow_multipleExamples(self.id(),
                                     lambda regionNumber=0, tempDir=None : getInputs(inputDir, mammalSequences),
                                     batchSystem=self.batchSystem,
                                     buildToilStats=True)

if __name__ == '__main__':
    unittest.main()
