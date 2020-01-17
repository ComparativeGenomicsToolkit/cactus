#!/usr/bin/env python3

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import os
import xml.etree.ElementTree as ET

from sonLib.bioio import TestStatus
from sonLib.bioio import getLogLevelString

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import runWorkflow_multipleExamples
from cactus.shared.test import silentOnSuccess

from cactus.shared.common import cactusRootPath
from sonLib.bioio import getTempFile

from cactus.shared.common import cactus_call

class TestCase(unittest.TestCase):
    @unittest.skip('too slow')
    @silentOnSuccess
    def testCactus_Random_Greedy(self):
        runCactus_Random(self, "greedy")

    @unittest.skip('too slow')
    @silentOnSuccess
    def testCactus_Random_Blossum(self):
        runCactus_Random(self, "blossom5")

    @unittest.skip('too slow')
    @silentOnSuccess
    def testCactus_Random_MaxCardinality(self):
        runCactus_Random(self, "maxCardinality")

    @unittest.skip('too slow')
    @silentOnSuccess
    def testCactus_Random_MaxWeight(self):
        runCactus_Random(self, "maxWeight")

    @unittest.skip('too slow')
    @silentOnSuccess
    def testCactus_Blanchette_Blossum(self):
        runCactus_Blanchette(self, "blossom5")

    def testCuTest(self):
        cactus_call(parameters=["referenceTests", getLogLevelString()])

def runCactus_Blanchette(self, matchingAlgorithm):
    configFile = getConfigFile(matchingAlgorithm)
    runWorkflow_multipleExamples(getCactusInputs_blanchette,
                                 testRestrictions=(TestStatus.TEST_SHORT,), inverseTestRestrictions=True,
                                 buildReference=True,
                                 configFile=configFile)
    os.remove(configFile)

def runCactus_Random(self, matchingAlgorithm):
    configFile = getConfigFile(matchingAlgorithm)
    runWorkflow_multipleExamples(getCactusInputs_random,
                                 testNumber=TestStatus.getTestSetup(),
                                 buildReference=True,
                                 configFile=configFile)
    os.remove(configFile)

def getConfigFile(matchingAlgorithm="greedy"):
    tempConfigFile = getTempFile(rootDir="./", suffix=".xml")
    config = ET.parse(os.path.join(cactusRootPath(), "cactus_progressive_config.xml")).getroot()
    #Set the matching algorithm
    config.find("reference").attrib["matching_algorithm"] = matchingAlgorithm
    #Now print the file..
    fileHandle = open(tempConfigFile, 'w')
    ET.ElementTree(config).write(fileHandle)
    fileHandle.close()
    return os.path.abspath(tempConfigFile)

if __name__ == '__main__':
    unittest.main()
