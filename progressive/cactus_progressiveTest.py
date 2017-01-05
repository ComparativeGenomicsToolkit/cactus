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

from sonLib.bioio import TestStatus
from sonLib.bioio import getTempFile
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import getRandomSequence
from sonLib.nxnewick import NXNewick

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import getCactusInputs_encode
from cactus.shared.test import getCactusInputs_chromosomeX
from cactus.shared.test import runWorkflow_multipleExamples
from cactus.shared.test import getBatchSystem
from cactus.shared.test import silentOnSuccess

from cactus.shared.experimentWrapper import ExperimentWrapper

from cactus.shared.common import cactusRootPath
from cactus.shared.common import runCactusProgressive
from cactus.shared.common import runCactusCreateMultiCactusProject
from cactus.shared.configWrapper import ConfigWrapper
from jobTree.src.common import runJobTreeStatusAndFailIfNotComplete

class TestCase(unittest.TestCase):
    def setUp(self):
        self.batchSystem = "singleMachine"
        if getBatchSystem() != None:
            self.batchSystem = getBatchSystem()
        unittest.TestCase.setUp(self)
        #Load the config file, turn on the checks.
        configWrapper = ConfigWrapper(ET.parse(os.path.join(cactusRootPath(), "cactus_progressive_config.xml")).getroot())
        configWrapper.turnAllModesOn()
        self.tempDir = getTempDirectory(os.getcwd())
        self.configFile = os.path.join(self.tempDir, "tempConfig.xml")
        configWrapper.writeXML(self.configFile)

    def tearDown(self):
        system("rm -rf %s" % self.tempDir)

    @silentOnSuccess
    def testCactus_Random(self):
        runWorkflow_multipleExamples(getCactusInputs_random,
                                     testNumber=2,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)

    @silentOnSuccess
    def testCactus_Random_UseSubtreeRoot(self):
        """Tests that cactus doesn't crash when aligning a subtree of a larger
        species tree."""
        runWorkflow_multipleExamples(getCactusInputs_random,
                                     testNumber=2,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveWithSubtreeRootFunction)

    @silentOnSuccess
    def testCactus_Blanchette(self):
        runWorkflow_multipleExamples(getCactusInputs_blanchette,
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_MEDIUM,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)

    @silentOnSuccess
    def testCactus_Encode(self):
        runWorkflow_multipleExamples(getCactusInputs_encode,
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_LONG,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)

    @silentOnSuccess
    def testCactus_Chromosomes(self):
        runWorkflow_multipleExamples(getCactusInputs_chromosomeX,
                                     testRestrictions=(TestStatus.TEST_VERY_LONG,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)

    def progressiveWithSubtreeRootFunction(self, experimentFile, jobTreeDir,
                          batchSystem, buildAvgs,
                          buildReference,
                          buildHal,
                          buildFasta,
                          jobTreeStats):
        """Choose an arbitrary subtree from the larger species tree to run the
        alignment on. This function is necessary to keep
        runWorkflow_multipleExamples general (specifying a subtree
        root doesn't make sense for runCactusWorkflow).
        """
        # Get valid internal nodes that are the root of the subtree we
        # want to align
        expWrapper = ExperimentWrapper(ET.parse(experimentFile).getroot())
        tree = expWrapper.getTree()
        validNodes = []
        for node in tree.postOrderTraversal():
            if tree.hasName(node) and not tree.isLeaf(node):
                validNodes.append(tree.getName(node))

        # Choose a random valid subtree root (NB: the entire species
        # tree is a valid subtree)
        subtreeRoot = random.choice(validNodes)
        logger.info("Chose subtree root %s to test from species tree "
                    "%s" % (subtreeRoot, NXNewick().writeString(tree)))

        self.progressiveFunction(experimentFile, jobTreeDir,
                                 batchSystem, buildAvgs,
                                 buildReference,
                                 buildHal,
                                 buildFasta,
                                 jobTreeStats, subtreeRoot)

    def progressiveFunction(self, experimentFile, jobTreeDir,
                            batchSystem, buildAvgs,
                            buildReference,
                            buildHal,
                            buildFasta,
                            jobTreeStats,
                            subtreeRoot=None):
        tempDir = getTempDirectory(os.getcwd())
        tempExperimentDir = os.path.join(tempDir, "exp")
        runCactusCreateMultiCactusProject(experimentFile,
                                          tempExperimentDir,
                                          fixNames=False,
                                          root=subtreeRoot)
        logger.info("Put the temporary files in %s" % tempExperimentDir)
        runCactusProgressive(os.path.join(tempExperimentDir, "exp_project.xml"),
                             jobTreeDir,
                             batchSystem=batchSystem,
                             buildAvgs=buildAvgs,
                             jobTreeStats=jobTreeStats)
        runJobTreeStatusAndFailIfNotComplete(jobTreeDir)
        system("rm -rf %s" % tempDir)
        
if __name__ == '__main__':
    unittest.main()
