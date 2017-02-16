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
import glob
import xml.etree.ElementTree as ET

from operator import itemgetter

from sonLib.bioio import TestStatus
from sonLib.bioio import getTempFile
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import getRandomSequence
from sonLib.bioio import fastaRead, fastaWrite
from sonLib.nxnewick import NXNewick

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import getCactusInputs_funkyHeaderNames
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
        #Load the config file, turn on the checks.
        configWrapper = ConfigWrapper(ET.parse(os.path.join(cactusRootPath(), "cactus_progressive_config.xml")).getroot())
        configWrapper.turnAllModesOn()
        configWrapper.turnOffHeaderChecks()
        self.tempDir = getTempDirectory(os.getcwd())
        self.configFile = os.path.join(self.tempDir, "tempConfig.xml")
        configWrapper.writeXML(self.configFile)

    def tearDown(self):
        system("rm -rf %s" % self.tempDir)

    @silentOnSuccess
    def testCactus_Random(self):
        # TODO: this doesn't actually need the test data, but this is
        # being used as a signal that the tester doesn't want to run
        # the long tests. The tests should be refactored soon.
        if "SON_TRACE_DATASETS" not in os.environ:
            return
        runWorkflow_multipleExamples(getCactusInputs_random,
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)

    @silentOnSuccess
    def testCactus_Random_UseSubtreeRoot(self):
        # TODO: this doesn't actually need the test data, but this is
        # being used as a signal that the tester doesn't want to run
        # the long tests. The tests should be refactored soon.
        if "SON_TRACE_DATASETS" not in os.environ:
            return
        """Tests that cactus doesn't crash when aligning a subtree of a larger
        species tree."""
        runWorkflow_multipleExamples(getCactusInputs_random,
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveWithSubtreeRootFunction)

    @silentOnSuccess
    def testCactus_ensureFunkyHeaderNamesArentMangled(self):
        """Ensure header names with characters like "|", " " aren't mangled."""
        if "SON_TRACE_DATASETS" not in os.environ:
            return
        runWorkflow_multipleExamples(getCactusInputs_funkyHeaderNames,
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)

    @silentOnSuccess
    def testCactus_Blanchette(self):
        if "SON_TRACE_DATASETS" not in os.environ:
            return
        runWorkflow_multipleExamples(getCactusInputs_blanchette,
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_MEDIUM,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)

    @silentOnSuccess
    def testCactus_Encode(self):
        if "SON_TRACE_DATASETS" not in os.environ:
            return
        runWorkflow_multipleExamples(getCactusInputs_encode,
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_LONG,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)

    @silentOnSuccess
    def testCactus_Chromosomes(self):
        if "SON_TRACE_DATASETS" not in os.environ:
            return
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

        # Check that the headers and sequences in the output are the
        # same as the sequences in the input (minus differences in
        # repeat-masking)
        exp = ExperimentWrapper(ET.parse(experimentFile).getroot())
        seqMap = exp.buildSequenceMap()
        # Maps genome name -> headers in fasta
        headers = {}
        for genomeName, inputSequencePath in seqMap.items():
            if os.path.isdir(inputSequencePath):
                # Some "input sequence paths" are actually provided as
                # directories containing multiple FASTAs
                concatenatedPath = getTempFile()
                system("cat %s/* > %s" % (inputSequencePath, concatenatedPath))
                inputSequencePath = concatenatedPath
            headers[genomeName] = list(map(itemgetter(0), fastaRead(inputSequencePath)))

        # check headers inside .c2h output
        for expPath in glob.glob('%s/*/*_experiment.xml' % (tempExperimentDir)):
            subExp = ExperimentWrapper(ET.parse(expPath).getroot())
            outgroups = subExp.getOutgroupEvents()
            c2hPath = subExp.getHALPath()
            with open(c2hPath) as f:
                for line in f:
                    fields = line.split('\t')
                    if fields[0] == 's':
                        # Sequence line
                        genome = fields[1][1:-1]
                        header = fields[2][1:-1]
                        if genome in headers and genome not in outgroups:
                            # This genome is an input genome
                            self.assertTrue(header in headers[genome],
                                            'Header %s from output c2h %s not found in input fa %s'
                                            ' for genome %s' % (header, c2hPath, seqMap[genome], genome))

        runJobTreeStatusAndFailIfNotComplete(jobTreeDir)
        system("rm -rf %s" % tempDir)
        
if __name__ == '__main__':
    unittest.main()
