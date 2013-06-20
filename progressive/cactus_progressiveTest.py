#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
"""Tests the core pipeline.
"""

import unittest
import os
import sys
import random

from cactus.shared.test import parseCactusSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger
from sonLib.bioio import system

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import getCactusInputs_encode
from cactus.shared.test import getCactusInputs_chromosomeX

from cactus.shared.test import runWorkflow_multipleExamples

from cactus.shared.test import getBatchSystem

from cactus.shared.common import runCactusProgressive
from cactus.shared.common import runCactusCreateMultiCactusProject
from jobTree.src.common import runJobTreeStatusAndFailIfNotComplete

class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.batchSystem = "singleMachine"
        if getBatchSystem() != None:
            self.batchSystem = getBatchSystem()
        unittest.TestCase.setUp(self)
        self.useOutgroup = False
        self.doSelfAlignment = False
        self.configFile = "defaultProgressive"
        if random.random() > 0.5:
            self.configFile = "defaultProgressiveFast"
        
    def testCactus_Random(self):
        runWorkflow_multipleExamples(getCactusInputs_random, 
                                     testNumber=2,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)
        
    def testCactus_Random_UseOutgroup(self):
        self.useOutgroup = True
        runWorkflow_multipleExamples(getCactusInputs_random, 
                                     testNumber=2,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)
        
    def testCactus_Random_DoSelfAlignment(self):
        self.doSelfAlignment = True
        runWorkflow_multipleExamples(getCactusInputs_random, 
                                     testNumber=2,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)
        
    def testCactus_Random_UseOutgroupAndDoSelfAlignment(self):
        self.useOutgroup = True
        self.doSelfAlignment = True
        runWorkflow_multipleExamples(getCactusInputs_random, 
                                     testNumber=2,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)
        
    def testCactus_Blanchette(self):
        runWorkflow_multipleExamples(getCactusInputs_blanchette, 
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_MEDIUM,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)
    
    def testCactus_Blanchette_UseOutgroupAndDoSelfAlignment(self):
        self.useOutgroup = True
        self.doSelfAlignment = True
        runWorkflow_multipleExamples(getCactusInputs_blanchette, 
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_MEDIUM,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)
                
    def testCactus_Encode(self): 
        runWorkflow_multipleExamples(getCactusInputs_encode, 
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_LONG,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)
    
    def testCactus_Chromosomes(self):
        runWorkflow_multipleExamples(getCactusInputs_chromosomeX, 
                                     testRestrictions=(TestStatus.TEST_VERY_LONG,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)
    
    def progressiveFunction(self, experimentFile, jobTreeDir, 
                          batchSystem, buildAvgs, 
                          buildReference,
                          buildHal,
                          buildMaf,
                          buildFasta,
                          jobTreeStats):
        tempDir = getTempDirectory(os.getcwd())
        tempExperimentDir = os.path.join(tempDir, "exp")
        runCactusCreateMultiCactusProject(experimentFile, 
                                          tempExperimentDir,
                                          fixNames = False)
        logger.info("Put the temporary files in %s" % tempExperimentDir)
        runCactusProgressive(os.path.join(tempExperimentDir, "exp_project.xml"), 
                             jobTreeDir, 
                             batchSystem=batchSystem, 
                             buildAvgs=buildAvgs,
                             jobTreeStats=jobTreeStats)
        runJobTreeStatusAndFailIfNotComplete(jobTreeDir)
        system("rm -rf %s" % tempDir)
    
def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
