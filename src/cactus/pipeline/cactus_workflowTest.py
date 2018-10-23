#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
"""Tests the core pipeline.
"""

import unittest
import os
import shutil
import xml.etree.ElementTree as ET
from tempfile import mkdtemp, NamedTemporaryFile
from textwrap import dedent

from sonLib.bioio import TestStatus, newickTreeParser, getTempFile

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_randomWithConstraints
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import getCactusInputs_encode
from cactus.shared.test import getCactusInputs_chromosomeX
from cactus.shared.test import runWorkflow_multipleExamples
from cactus.shared.test import getBatchSystem
from cactus.shared.test import silentOnSuccess
from cactus.shared.test import initialiseGlobalDatabaseConf

from cactus.shared.common import cactusRootPath

from cactus.pipeline.cactus_workflow import getOptionalAttrib, extractNode, findRequiredNode, \
    getJobNode, CactusJob, getLongestPath, inverseJukesCantor, \
    CactusSetReferenceCoordinatesDownRecursion, prependUniqueIDs

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
    @unittest.skip("")
    def testCactus_random(self):
        runWorkflow_multipleExamples(getCactusInputs_random, 
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     buildAvgs=True, buildReference=True,
                                     batchSystem=self.batchSystem, buildToilStats=True)

    @silentOnSuccess
    @unittest.skip("")
    def testCactus_randomWithConstraints(self):
        runWorkflow_multipleExamples(getCactusInputs_randomWithConstraints, 
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_SHORT,),
                                     buildAvgs=True, buildReference=True,
                                     batchSystem=self.batchSystem, buildToilStats=True,
                                     useConstraints=True)

    @silentOnSuccess
    @unittest.skip("")
    def testCactus_blanchette(self):
        runWorkflow_multipleExamples(getCactusInputs_blanchette, 
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_MEDIUM,),
                                     buildAvgs=True, buildReference=True,
                                     batchSystem=self.batchSystem, buildToilStats=True)
    @silentOnSuccess
    @unittest.skip("")
    def testCactus_encode(self): 
        runWorkflow_multipleExamples(getCactusInputs_encode, 
                                     testNumber=1,
                                     testRestrictions=(TestStatus.TEST_LONG,),
                                     buildAvgs=True, buildReference=True,
                                     batchSystem=self.batchSystem, buildToilStats=True)
    @silentOnSuccess
    @unittest.skip("")
    def testCactus_chromosomes(self):
        runWorkflow_multipleExamples(getCactusInputs_chromosomeX, 
                                     testRestrictions=(TestStatus.TEST_VERY_LONG,),
                                     batchSystem=self.batchSystem, buildToilStats=True)

    @silentOnSuccess
    def testCactus_splitBarJobs(self):
        """Exercise the code paths in bar that only occur on large jobs."""
        # Modify the bar node in the config file so that
        # cactus_workflow will split bar jobs even on this small
        # example
        initialiseGlobalDatabaseConf('<st_kv_database_conf type="kyoto_tycoon"><kyoto_tycoon in_memory="1" port="1978" snapshot="0"/></st_kv_database_conf>')
        tempConfigFile = getTempFile()
        tempConfigTree = ET.parse(self.configFile)
        tempConfigNode = tempConfigTree.getroot()
        tempConfigNode.find("bar").find("CactusBarWrapper").set("maxFlowerGroupSize", "10")
        tempConfigNode.find("bar").find("CactusBarWrapperLarge").set("maxFlowerGroupSize", "10")
        tempConfigNode.find("bar").set("veryLargeEndSize", "20")
        tempConfigNode.find("bar").set("largeEndSize", "10")
        tempConfigNode.find("bar").set("bandingLimit", "5")
        tempConfigTree.write(tempConfigFile)
        runWorkflow_multipleExamples(getCactusInputs_random,
                                     testNumber=1,
                                     batchSystem=self.batchSystem,
                                     configFile=tempConfigFile)
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
        self.assertEquals(findRequiredNode(self.configNode, "caf"), self.configNode.findall("caf")[0])

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
        self.assertEquals(job.getOptionalPhaseAttrib("diagonalExpansion", typeFn=int), 20)
        self.assertEquals(job.getOptionalPhaseAttrib("doesntExist", typeFn=int, default=1), 1)
        self.assertEquals(job.getOptionalJobAttrib("memory", typeFn=int), 10)
        self.assertEquals(job.getOptionalJobAttrib("cpu", typeFn=int, default=1), 2)
        self.assertEquals(job.getOptionalJobAttrib("overlargeCpu", typeFn=int, default=-1), -1)

    def testGetLongestPath(self):
        self.assertAlmostEquals(getLongestPath(newickTreeParser("(b(a:0.5):0.5,b(a:1.5):0.5)")), 2.0)
        self.assertAlmostEquals(getLongestPath(newickTreeParser("(b(a:0.5):0.5,b(a:1.5,c:10):0.5)")), 10.5)
        self.assertAlmostEquals(getLongestPath(newickTreeParser("(b(a:0.5):0.5,b(a:1.5,c:10,e,f:20):0.5)")), 20.5)

    def testInverseJukesCantor(self):
        self.assertAlmostEquals(inverseJukesCantor(0.5), 0.36493716072555599)
        self.assertAlmostEquals(inverseJukesCantor(1.0), 0.55230214641320496)
        self.assertAlmostEquals(inverseJukesCantor(10.0), 0.74999878530240571)
        self.assertAlmostEquals(inverseJukesCantor(100000.0), 0.75)

    def testPrependUniqueIDs(self):
        # Create fake FASTA files with some interesting headers.
        with NamedTemporaryFile() as fasta1, NamedTemporaryFile() as fasta2:
            fasta1.write(dedent("""
            >C10856240  2.0
            ACTAGAGG
            G

            GG
            >foo bar  baz
            ACTGACGATgacgat
            >emptyseq
            """))
            fasta2.write(dedent("""
            > space
            GTGC
            >id=1||header
            ATCC
            """))
            fasta1.flush()
            fasta2.flush()
            outDir = mkdtemp()
            outputPaths = prependUniqueIDs([fasta1.name, fasta2.name], outDir)

        assert len(outputPaths) == 2
        with open(outputPaths[0]) as f:
            self.assertEquals(f.read(), dedent("""
            >id=0|C10856240  2.0
            ACTAGAGG
            G

            GG
            >id=0|foo bar  baz
            ACTGACGATgacgat
            >id=0|emptyseq
            """))
        with open(outputPaths[1]) as f:
            self.assertEquals(f.read(), dedent("""
            >id=1| space
            GTGC
            >id=1|id=1||header
            ATCC
            """))
        shutil.rmtree(outDir)

if __name__ == '__main__':
    unittest.main()
