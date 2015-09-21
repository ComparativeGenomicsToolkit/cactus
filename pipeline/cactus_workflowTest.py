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
from sonLib.bioio import TestStatus, newickTreeParser

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_randomWithConstraints
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import getCactusInputs_encode
from cactus.shared.test import getCactusInputs_chromosomeX

from cactus.shared.test import runWorkflow_multipleExamples

from cactus.shared.test import getBatchSystem

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
        
    def testCactus_random(self):
        runWorkflow_multipleExamples(getCactusInputs_random, 
                                     testNumber=5,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     buildAvgs=True, buildReference=True,
                                     batchSystem=self.batchSystem, buildToilStats=True)
        
    def testCactus_randomWithConstraints(self):
        runWorkflow_multipleExamples(getCactusInputs_randomWithConstraints, 
                                     testNumber=5,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     buildAvgs=True, buildReference=True,
                                     batchSystem=self.batchSystem, buildToilStats=True,
                                     useConstraints=True)
        
    def testCactus_blanchette(self):
        runWorkflow_multipleExamples(getCactusInputs_blanchette, 
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_MEDIUM,),
                                     buildAvgs=True, buildReference=True,
                                     batchSystem=self.batchSystem, buildToilStats=True)
                
    def testCactus_encode(self): 
        runWorkflow_multipleExamples(getCactusInputs_encode, 
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_LONG,),
                                     buildAvgs=True, buildReference=True,
                                     batchSystem=self.batchSystem, buildToilStats=True)
    
    def testCactus_chromosomes(self):
        runWorkflow_multipleExamples(getCactusInputs_chromosomeX, 
                                     testRestrictions=(TestStatus.TEST_VERY_LONG,),
                                     batchSystem=self.batchSystem, buildToilStats=True)
        
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
        
    def testGetJobNode(self):
        class CactusTestJob(CactusJob):
            pass
        class CactusTestJob2(CactusJob):
            pass
        node = ET.SubElement(self.barNode, "CactusTestJob")
        self.assertEquals(node, getJobNode(self.barNode, CactusTestJob))
        self.assertEquals(None, getJobNode(self.barNode, CactusTestJob2))
        node2 = ET.SubElement(self.barNode, "CactusSetReferenceCoordinatesDownRecursion")
        self.assertEquals(node2, getJobNode(self.barNode, CactusSetReferenceCoordinatesDownRecursion))
    
    def testCactusJob(self):
        class CactusTestJob(CactusJob):
            pass
        node = ET.SubElement(self.barNode, "CactusTestJob", attrib={ "memory":10, "cpu":2,  "overlargeMemory":20 })
        job = CactusTestJob(self.barNode, self.barNode)
        self.assertEquals(job.jobNode, node)
        self.assertEquals(job.memory, 10)
        self.assertEquals(job.cores, 2)
        job = CactusTestJob(self.barNode, self.barNode, overlarge=True)
        self.assertEquals(job.memory, 20)
        #self.assertEquals(job.cores, sys.maxint)
        self.assertEquals(job.getOptionalPhaseAttrib("diagonalExpansion", typeFn=int), 20)
        self.assertEquals(job.getOptionalPhaseAttrib("doesntExist", typeFn=int, default=1), 1)
        self.assertEquals(job.getOptionalJobAttrib("memory", typeFn=int), 10)
        self.assertEquals(job.getOptionalJobAttrib("cpu", typeFn=int, default=1), 2)
        self.assertEquals(job.getOptionalJobAttrib("overlargeCpu", typeFn=int, default=-1), -1)
        class CactusTestJob2(CactusJob):
            pass
        job = CactusTestJob2(self.barNode, self.barNode)
        self.assertEquals(job.jobNode, None)
        #self.assertEquals(job.memory, sys.maxint)
        #self.assertEquals(job.cores, sys.maxint)
    
    def testGetLongestPath(self):
        self.assertAlmostEquals(getLongestPath(newickTreeParser("(b(a:0.5):0.5,b(a:1.5):0.5)")), 2.0)
        self.assertAlmostEquals(getLongestPath(newickTreeParser("(b(a:0.5):0.5,b(a:1.5,c:10):0.5)")), 10.5)
        self.assertAlmostEquals(getLongestPath(newickTreeParser("(b(a:0.5):0.5,b(a:1.5,c:10,e,f:20):0.5)")), 20.5)
    
    def testInverseJukesCantor(self):
        self.assertAlmostEquals(inverseJukesCantor(0.5), 0.36493716072555599)
        self.assertAlmostEquals(inverseJukesCantor(1.0), 0.55230214641320496)
        self.assertAlmostEquals(inverseJukesCantor(10.0), 0.74999878530240571)
        self.assertAlmostEquals(inverseJukesCantor(100000.0), 0.75)

def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    from cactus.pipeline.cactus_workflowTest import *
    main()
