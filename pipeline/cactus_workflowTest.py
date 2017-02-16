#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
"""Tests the core pipeline.
"""

import unittest
import os
import sys

from sonLib.bioio import TestStatus, newickTreeParser, getTempFile

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_randomWithConstraints
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import getCactusInputs_encode
from cactus.shared.test import getCactusInputs_chromosomeX
from cactus.shared.test import runWorkflow_multipleExamples
from cactus.shared.test import getBatchSystem
from cactus.shared.test import silentOnSuccess

from cactus.shared.common import cactusRootPath

from cactus.pipeline.cactus_workflow import *

class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.batchSystem = "singleMachine"
        if getBatchSystem() != None:
            self.batchSystem = getBatchSystem()
        unittest.TestCase.setUp(self)
        self.configFile = os.path.join(cactusRootPath(), "cactus_config.xml")
        self.configNode = ET.parse(self.configFile).getroot()
        self.barNode = self.configNode.find("bar")
        assert self.barNode != None

    @silentOnSuccess
    def testCactus_random(self):
        runWorkflow_multipleExamples(getCactusInputs_random, 
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     buildAvgs=True, buildReference=True,
                                     batchSystem=self.batchSystem, buildJobTreeStats=True)

    @silentOnSuccess
    def testCactus_randomWithConstraints(self):
        runWorkflow_multipleExamples(getCactusInputs_randomWithConstraints, 
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     buildAvgs=True, buildReference=True,
                                     batchSystem=self.batchSystem, buildJobTreeStats=True,
                                     useConstraints=True)

    @silentOnSuccess
    def testCactus_blanchette(self):
        runWorkflow_multipleExamples(getCactusInputs_blanchette, 
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_MEDIUM,),
                                     buildAvgs=True, buildReference=True,
                                     batchSystem=self.batchSystem, buildJobTreeStats=True)

    @silentOnSuccess
    def testCactus_encode(self): 
        runWorkflow_multipleExamples(getCactusInputs_encode, 
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_LONG,),
                                     buildAvgs=True, buildReference=True,
                                     batchSystem=self.batchSystem, buildJobTreeStats=True)

    @silentOnSuccess
    def testCactus_chromosomes(self):
        runWorkflow_multipleExamples(getCactusInputs_chromosomeX, 
                                     testRestrictions=(TestStatus.TEST_VERY_LONG,),
                                     batchSystem=self.batchSystem, buildJobTreeStats=True)
    @silentOnSuccess
    def testCactus_splitBarJobs(self):
        """Exercise the code paths in bar that only occur on large jobs."""
        # Modify the bar node in the config file so that
        # cactus_workflow will split bar jobs even on this small
        # example
        tempConfigFile = getTempFile()
        tempConfigTree = ET.parse(self.configFile)
        tempConfigNode = tempConfigTree.getroot()
        tempConfigNode.find("bar").find("CactusBarWrapper").set("maxFlowerGroupSize", "10")
        tempConfigNode.find("bar").set("veryLargeEndSize", "0")
        tempConfigNode.find("bar").set("largeEndSize", "0")
        tempConfigTree.write(tempConfigFile)
        runWorkflow_multipleExamples(getCactusInputs_blanchette,
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_LONG,),
                                     batchSystem=self.batchSystem,
                                     configFile=tempConfigFile, buildJobTreeStats=True)
        os.remove(tempConfigFile)

    def testGetOptionalAttrib(self):
        self.assertEquals("0", getOptionalAttrib(self.barNode, "minimumBlockDegree"))
        self.assertEquals(0, getOptionalAttrib(self.barNode, "minimumBlockDegree", typeFn=int, default=1))
        self.assertEquals(None, getOptionalAttrib(self.barNode, "doesntExist"))
        self.assertEquals(1, getOptionalAttrib(self.barNode, "doesntExist", typeFn=int, default=1))
    
    def testFindRequiredNode(self):
        self.assertEquals(findRequiredNode(self.configNode, "bar"), self.barNode)
        try:
            findRequiredNode(self.configNode, "doesntExist")
            self.assertTrue(0)
        except:
            pass
        self.assertEquals(findRequiredNode(self.configNode, "caf", index=0), self.configNode.findall("caf")[0])
        try:
            findRequiredNode(self.configNode, "caf", index=1)
            self.assertTrue(0)
        except:
            pass
    
    def testExtractNode(self):
        subNode = ET.SubElement(self.barNode, "CactusSetReferenceCoordinatesDownRecursion", { "memory":"10" })
        barNodeCopy = extractNode(self.barNode)
        barNodeCopy.attrib["added"] = "1"
        self.assertFalse("added" in self.barNode.attrib)
        self.barNode.attrib["added2"] = "1"
        self.assertTrue("added2" in self.barNode.attrib)
        self.assertFalse("added2" in barNodeCopy.attrib)
        self.assertEquals(subNode, self.barNode.find("CactusSetReferenceCoordinatesDownRecursion"))
        subNodeCopy = barNodeCopy.find("CactusSetReferenceCoordinatesDownRecursion")
        self.assertTrue(subNodeCopy != None)
        self.assertEquals("10", subNodeCopy.attrib["memory"])
        
    def testGetTargetNode(self):
        class CactusTestTarget(CactusTarget):
            pass
        class CactusTestTarget2(CactusTarget):
            pass
        node = ET.SubElement(self.barNode, "CactusTestTarget")
        self.assertEquals(node, getTargetNode(self.barNode, CactusTestTarget))
        self.assertEquals(None, getTargetNode(self.barNode, CactusTestTarget2))
        node2 = ET.SubElement(self.barNode, "CactusSetReferenceCoordinatesDownRecursion")
        self.assertEquals(node2, getTargetNode(self.barNode, CactusSetReferenceCoordinatesDownRecursion))
    
    def testCactusTarget(self):
        class CactusTestTarget(CactusTarget):
            pass
        node = ET.SubElement(self.barNode, "CactusTestTarget", attrib={ "memory":10, "cpu":2,  "overlargeMemory":20 })
        target = CactusTestTarget(self.barNode, self.barNode)
        self.assertEquals(target.targetNode, node)
        self.assertEquals(target.getMemory(), 10)
        self.assertEquals(target.getCpu(), 2)
        target = CactusTestTarget(self.barNode, self.barNode, overlarge=True)
        self.assertEquals(target.getMemory(), 20)
        self.assertEquals(target.getCpu(), sys.maxint)
        self.assertEquals(target.getOptionalPhaseAttrib("diagonalExpansion", typeFn=int), 20)
        self.assertEquals(target.getOptionalPhaseAttrib("doesntExist", typeFn=int, default=1), 1)
        self.assertEquals(target.getOptionalTargetAttrib("memory", typeFn=int), 10)
        self.assertEquals(target.getOptionalTargetAttrib("cpu", typeFn=int, default=1), 2)
        self.assertEquals(target.getOptionalTargetAttrib("overlargeCpu", typeFn=int, default=-1), -1)
        class CactusTestTarget2(CactusTarget):
            pass
        target = CactusTestTarget2(self.barNode, self.barNode)
        self.assertEquals(target.targetNode, None)
        self.assertEquals(target.getMemory(), sys.maxint)
        self.assertEquals(target.getCpu(), sys.maxint)
    
    def testGetLongestPath(self):
        self.assertAlmostEquals(getLongestPath(newickTreeParser("(b(a:0.5):0.5,b(a:1.5):0.5)")), 2.0)
        self.assertAlmostEquals(getLongestPath(newickTreeParser("(b(a:0.5):0.5,b(a:1.5,c:10):0.5)")), 10.5)
        self.assertAlmostEquals(getLongestPath(newickTreeParser("(b(a:0.5):0.5,b(a:1.5,c:10,e,f:20):0.5)")), 20.5)

    def testInverseJukesCantor(self):
        self.assertAlmostEquals(inverseJukesCantor(0.5), 0.36493716072555599)
        self.assertAlmostEquals(inverseJukesCantor(1.0), 0.55230214641320496)
        self.assertAlmostEquals(inverseJukesCantor(10.0), 0.74999878530240571)
        self.assertAlmostEquals(inverseJukesCantor(100000.0), 0.75)

if __name__ == '__main__':
    from cactus.pipeline.cactus_workflowTest import *
    unittest.main()
