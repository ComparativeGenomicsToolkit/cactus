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
import xml.etree.ElementTree as ET

from cactus.shared.test import parseCactusSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempFile
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import getRandomSequence

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import getCactusInputs_encode
from cactus.shared.test import getCactusInputs_chromosomeX

from cactus.shared.test import runWorkflow_multipleExamples

from cactus.shared.test import getBatchSystem

from cactus.shared.common import cactusRootPath
from cactus.shared.common import runCactusProgressive
from cactus.shared.common import runCactusCreateMultiCactusProject
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import runToilStatusAndFailIfNotComplete

class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.batchSystem = "singleMachine"
        if getBatchSystem() != None:
            self.batchSystem = getBatchSystem()
        unittest.TestCase.setUp(self)
        self.useOutgroup = False
        self.doSelfAlignment = False
        #Load the config file, turn on the checks.
        configWrapper = ConfigWrapper(ET.parse(os.path.join(cactusRootPath(), "cactus_progressive_config.xml")).getroot())
        configWrapper.turnAllModesOn()
        self.tempDir = getTempDirectory(os.getcwd())
        self.configFile = os.path.join(self.tempDir, "tempConfig.xml")
        configWrapper.writeXML(self.configFile)
    
    def tearDown(self):
        system("rm -rf %s" % self.tempDir)
        
    def testCactus_Random(self):
        runWorkflow_multipleExamples(getCactusInputs_random, 
                                     testNumber=2,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     batchSystem=self.batchSystem, buildToilStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)
        
    def testCactus_Random_UseOutgroup(self):
        self.useOutgroup = True
        runWorkflow_multipleExamples(getCactusInputs_random, 
                                     testNumber=2,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     batchSystem=self.batchSystem, buildToilStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)

    def testCactus_Random_UseRootOutgroup(self):
        runWorkflow_multipleExamples(getCactusInputs_random, 
                                     testNumber=2,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     batchSystem=self.batchSystem, buildToilStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveWithRootOutgroupFunction)
        
    def testCactus_Random_DoSelfAlignment(self):
        self.doSelfAlignment = True
        runWorkflow_multipleExamples(getCactusInputs_random, 
                                     testNumber=2,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     batchSystem=self.batchSystem, buildToilStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)
        
    def testCactus_Random_UseOutgroupAndDoSelfAlignment(self):
        self.useOutgroup = True
        self.doSelfAlignment = True
        runWorkflow_multipleExamples(getCactusInputs_random, 
                                     testNumber=2,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     batchSystem=self.batchSystem, buildToilStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)
        
    def testCactus_Blanchette(self):
        runWorkflow_multipleExamples(getCactusInputs_blanchette, 
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_MEDIUM,),
                                     batchSystem=self.batchSystem, buildToilStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)
    
    def testCactus_Blanchette_UseOutgroupAndDoSelfAlignment(self):
        self.useOutgroup = True
        self.doSelfAlignment = True
        runWorkflow_multipleExamples(getCactusInputs_blanchette, 
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_MEDIUM,),
                                     batchSystem=self.batchSystem, buildToilStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)
                
    def testCactus_Encode(self): 
        runWorkflow_multipleExamples(getCactusInputs_encode, 
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_LONG,),
                                     batchSystem=self.batchSystem, buildToilStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)
    
    def testCactus_Chromosomes(self):
        runWorkflow_multipleExamples(getCactusInputs_chromosomeX, 
                                     testRestrictions=(TestStatus.TEST_VERY_LONG,),
                                     batchSystem=self.batchSystem, buildToilStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)
    
    def progressiveWithRootOutgroupFunction(self, experimentFile, toilDir, 
                          batchSystem, buildAvgs, 
                          buildReference,
                          buildHal, 
                          buildFasta,
                          toilStats):
        """Add in a (random, small) root outgroup before calling
        progressiveFunction. This function is necessary to keep
        runWorkflow_multipleExamples general (root outgroups don't
        make sense for runCactusWorkflow).
        """
        rootOutgroupPath = getTempFile()
        outgroupSequence = getRandomSequence(length=random.choice(xrange(1, 10000)))
        self.progressiveFunction(experimentFile, toilDir,
                                 batchSystem, buildAvgs, 
                                 buildReference,
                                 buildHal, 
                                 buildFasta,
                                 toilStats, rootOutgroupPath, 1.0)
        system("rm -f %s" % rootOutgroupPath)

    def progressiveFunction(self, experimentFile, toilDir, 
                          batchSystem, buildAvgs, 
                          buildReference,
                          buildHal, 
                          buildFasta,
                          toilStats, rootOutgroupPath=None,
                          rootOutgroupDist=None):
        tempDir = getTempDirectory(os.getcwd())
        tempExperimentDir = os.path.join(tempDir, "exp")
        runCactusCreateMultiCactusProject(experimentFile, 
                                          tempExperimentDir,
                                          fixNames = False,
                                          rootOutgroupPath=rootOutgroupPath,
                                          rootOutgroupDist=rootOutgroupDist)
        logger.info("Put the temporary files in %s" % tempExperimentDir)
        runCactusProgressive(os.path.join(tempExperimentDir, "exp_project.xml"), 
                             toilDir, 
                             batchSystem=batchSystem, 
                             buildAvgs=buildAvgs,
                             toilStats=toilStats)
        runToilStatusAndFailIfNotComplete(toilDir)
        system("rm -rf %s" % tempDir)
    
def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
