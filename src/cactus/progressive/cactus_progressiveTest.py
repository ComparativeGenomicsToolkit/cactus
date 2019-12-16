#!/usr/bin/env python33

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
"""Tests the core pipeline.
"""

import unittest
import pytest
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
from cactus.shared.test import getCactusWorkflowExperimentForTest

from cactus.shared.experimentWrapper import ExperimentWrapper

from cactus.shared.common import cactusRootPath
from cactus.shared.common import runCactusProgressive
from cactus.progressive.cactus_createMultiCactusProject import runCreateMultiCactusProject
from cactus.shared.configWrapper import ConfigWrapper

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
        configWrapper.turnOffHeaderChecks()
        self.tempDir = getTempDirectory(os.getcwd())
        self.configFile = os.path.join(self.tempDir, "tempConfig.xml")
        configWrapper.writeXML(self.configFile)

    def tearDown(self):
        system("rm -rf %s" % self.tempDir)

    @silentOnSuccess
    @unittest.skip("")
    def testCactus_Random(self):
        # TODO: this doesn't actually need the test data, but this is
        # being used as a signal that the tester doesn't want to run
        # the long tests. The tests should be refactored soon.
        if "SON_TRACE_DATASETS" not in os.environ:
            return
        runWorkflow_multipleExamples(getCactusInputs_random,
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     batchSystem=self.batchSystem, buildToilStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)

    @silentOnSuccess
    @unittest.skip("")
    def testCactus_Random_UseSubtreeRoot(self):
        # TODO: this doesn't actually need the test data, but this is
        # being used as a signal that the tester doesn't want to run
        # the long tests. The tests should be refactored soon.
        if "SON_TRACE_DATASETS" not in os.environ:
            return
        """Tests that cactus doesn't crash when aligning a subtree of a larger
        species tree."""
        getBigSpeciesTree = lambda regionNumber=0,tempDir=None: getCactusInputs_random(regionNumber, tempDir, treeLeafNumber=5)
        runWorkflow_multipleExamples(getBigSpeciesTree,
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     batchSystem=self.batchSystem, buildToilStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveWithSubtreeRootFunction)

    @pytest.mark.skipif(os.environ.get("TRAVIS") is not None, reason="Travis doesn't have enough cores anymore")
    @silentOnSuccess
    def testCactus_Random_fixedAncestor(self):
        """Tests that cactus doesn't crash when aligning to a fixed ancestral sequence."""
        sequences, _ = getCactusInputs_random(treeLeafNumber=3)
        rootSeq = sequences.pop()
        # Create a star tree
        tree = '(%s)root;' % ",".join([str(x) + ":1.0" for x in range(len(sequences))])
        outputDir = getTempDirectory()
        experiment = getCactusWorkflowExperimentForTest(sequences, tree,
                                                        outputDir,
                                                        progressive=True)
        experiment.setSequenceID("root", rootSeq)
        experiment.setRootReconstructed(False)
        experimentFile = os.path.join(outputDir, "experiment.xml")
        experiment.writeXML(experimentFile)

        jobTreeDir = os.path.join(outputDir, "jobTree")

        self.progressiveFunction(experimentFile, jobTreeDir, 'singleMachine',
                                 False, True, True, False)

    @silentOnSuccess
    @unittest.skip("")
    def testCactus_ensureFunkyHeaderNamesArentMangled(self):
        """Ensure header names with characters like "|", " " aren't mangled."""
        if "SON_TRACE_DATASETS" not in os.environ:
            return
        runWorkflow_multipleExamples(getCactusInputs_funkyHeaderNames,
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     batchSystem=self.batchSystem, buildToilStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)

    @silentOnSuccess
    @unittest.skip("")
    def testCactus_Blanchette(self):
        runWorkflow_multipleExamples(getCactusInputs_blanchette,
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_MEDIUM,),
                                     batchSystem=self.batchSystem, buildToilStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)

    @silentOnSuccess
    @unittest.skip("")
    def testCactus_Encode(self):
        if "SON_TRACE_DATASETS" not in os.environ:
            return
        runWorkflow_multipleExamples(getCactusInputs_encode,
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_LONG,),
                                     batchSystem=self.batchSystem, buildToilStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)
    @silentOnSuccess
    @unittest.skip("")
    def testCactus_Chromosomes(self):
        if "SON_TRACE_DATASETS" not in os.environ:
            return
        runWorkflow_multipleExamples(getCactusInputs_chromosomeX,
                                     testRestrictions=(TestStatus.TEST_VERY_LONG,),
                                     batchSystem=self.batchSystem, buildToilStats=True,
                                     progressive=True,
                                     configFile=self.configFile,
                                     cactusWorkflowFunction=self.progressiveFunction)

    def progressiveWithSubtreeRootFunction(self, experimentFile, toilDir,
                                           batchSystem, buildAvgs,
                                           buildReference,
                                           buildHal,
                                           buildFasta,
                                           toilStats):
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
            if tree.hasName(node) and not tree.isLeaf(node) and tree.hasParent(node):
                validNodes.append(tree.getName(node))

        # Choose a random valid subtree root (excluding the species tree root)
        subtreeRoot = random.choice(validNodes)

        self.progressiveFunction(experimentFile, toilDir,
                                 batchSystem, buildAvgs,
                                 buildHal,
                                 buildFasta,
                                 toilStats, subtreeRoot)

    def progressiveFunction(self, experimentFile, toilDir,
                            batchSystem, buildAvgs,
                            buildHal,
                            buildFasta,
                            toilStats,
                            subtreeRoot=None):
        eW = ExperimentWrapper(ET.parse(experimentFile).getroot())
        seqFile = getTempFile()
        with open(seqFile, 'w') as f:
            tree = eW.getTree()
            newick = NXNewick().writeString(tree)
            f.write('%s\n' % newick)
            for genome in eW.getGenomesWithSequence():
                f.write('%s %s\n' % (genome, eW.getSequenceID(genome)))
        config = eW.getConfigPath()
        runCactusProgressive(seqFile,
                             config,
                             toilDir,
                             batchSystem=batchSystem,
                             buildAvgs=buildAvgs,
                             toilStats=toilStats)
        
if __name__ == '__main__':
    unittest.main()
