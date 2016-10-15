#!/usr/bin/env python

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
from cactus.shared.test import silentOnSuccess

class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.batchSystem = "singleMachine"
        if getBatchSystem() != None:
            self.batchSystem = getBatchSystem()
        unittest.TestCase.setUp(self)

    @silentOnSuccess
    def testEvolver_Primates_Loci1(self):
        if "SON_TRACE_DATASETS" not in os.environ:
            return
        inputDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "primates", "loci1")
        primateSequences = ("simChimp.chr6", "simGorilla.chr6", "simHuman.chr6", "simOrang.chr6")
        runWorkflow_multipleExamples(lambda regionNumber=0, tempDir=None : getInputs(inputDir, primateSequences),
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     batchSystem=self.batchSystem,
                                     buildJobTreeStats=True)

    @silentOnSuccess
    def testEvolver_Mammals_Loci1(self):
        if "SON_TRACE_DATASETS" not in os.environ:
            return
        inputDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "mammals", "loci1")
        mammalSequences = ("simCow.chr6", "simDog.chr6", "simHuman.chr6", "simMouse.chr6", "simRat.chr6")
        runWorkflow_multipleExamples(lambda regionNumber=0, tempDir=None : getInputs(inputDir, mammalSequences),
                                     testRestrictions=(TestStatus.TEST_MEDIUM,),
                                     batchSystem=self.batchSystem,
                                     buildJobTreeStats=True)

    @silentOnSuccess
    def testEvolver_Primates_Small(self):
        if "SON_TRACE_DATASETS" not in os.environ:
            return
        inputDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "primates", "small")
        primateSequences = ("simChimp.fa", "simGorilla.fa", "simHuman.fa", "simOrang.fa")
        runWorkflow_multipleExamples(lambda regionNumber=0, tempDir=None : getInputs(inputDir, primateSequences),
                                     testRestrictions=(TestStatus.TEST_MEDIUM,),
                                     batchSystem=self.batchSystem,
                                     buildJobTreeStats=True)

    @silentOnSuccess
    def testEvolver_Mammals_Medium(self):
        if "SON_TRACE_DATASETS" not in os.environ:
            return
        inputDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "mammals", "medium")
        mammalSequences = ("simCow.masked.fa", "simDog.masked.fa", "simHuman.masked.fa", "simMouse.masked.fa", "simRat.masked.fa")
        runWorkflow_multipleExamples(lambda regionNumber=0, tempDir=None : getInputs(inputDir, mammalSequences),
                                     testRestrictions=(TestStatus.TEST_LONG,),
                                     batchSystem=self.batchSystem,
                                     buildJobTreeStats=True)

    @silentOnSuccess
    def testEvolver_Primates_Large(self):
        if "SON_TRACE_DATASETS" not in os.environ:
            return
        inputDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "primates", "large")
        primateSequences = ("simChimp.masked.fa", "simGorilla.masked.fa", "simHuman.masked.fa", "simOrang.masked.fa")
        runWorkflow_multipleExamples(lambda regionNumber=0, tempDir=None : getInputs(inputDir, primateSequences),
                                     testRestrictions=(TestStatus.TEST_VERY_LONG,),
                                     batchSystem=self.batchSystem,
                                     buildJobTreeStats=True)

    @silentOnSuccess
    def testEvolver_Mammals_Large(self):
        if "SON_TRACE_DATASETS" not in os.environ:
            return
        inputDir = os.path.join(TestStatus.getPathToDataSets(), "evolver", "mammals", "large")
        mammalSequences = ("simCow.masked.fa", "simDog.masked.fa", "simHuman.masked.fa", "simMouse.masked.fa", "simRat.masked.fa")
        runWorkflow_multipleExamples(lambda regionNumber=0, tempDir=None : getInputs(inputDir, mammalSequences),
                                     testRestrictions=(TestStatus.TEST_VERY_LONG,),
                                     batchSystem=self.batchSystem,
                                     buildJobTreeStats=True)

if __name__ == '__main__':
    unittest.main()
