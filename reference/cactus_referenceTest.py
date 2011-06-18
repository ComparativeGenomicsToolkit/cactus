#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys
import os
import xml.etree.ElementTree as ET

from cactus.shared.test import parseCactusSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import system
from sonLib.bioio import getLogLevelString

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import runWorkflow_multipleExamples

from cactus.shared.common import cactusRootPath
from sonLib.bioio import getTempFile

class TestCase(unittest.TestCase):
    def testCactus_Random_Greedy(self):
        testCactus_Random(self, "greedy")
        
    def testCactus_Random_Blossum(self):
        testCactus_Random(self, "blossom5")
        
    def testCactus_Random_MaxCardinality(self):
        testCactus_Random(self, "maxCardinality")
    
    def testCactus_Random_MaxWeight(self):
        testCactus_Random(self, "maxWeight")
        
    def testCactus_Blanchette_Greedy(self):
        testCactus_Blanchette(self, "greedy")
        
    def testCactus_Blanchette_Blossum(self):
        testCactus_Blanchette(self, "blossom5")
        
    def testCactus_Blanchette_MaxCardinality(self):
        testCactus_Blanchette(self, "maxCardinality")
    
    def testCactus_Blanchette_MaxWeight(self):
        testCactus_Blanchette(self, "maxWeight")
        
    def testCuTest(self):
        system("referenceTests %s" % getLogLevelString())
            
def testCactus_Blanchette(self, matchingAlgorithm):
    configFile = getConfigFile(matchingAlgorithm)
    runWorkflow_multipleExamples(getCactusInputs_blanchette, 
                                 testRestrictions=(TestStatus.TEST_SHORT,), inverseTestRestrictions=True, 
                                 buildTrees=False, buildFaces=False, buildReference=True,
                                 configFile=configFile)
    os.remove(configFile)

def testCactus_Random(self, matchingAlgorithm):
    configFile = getConfigFile(matchingAlgorithm)
    runWorkflow_multipleExamples(getCactusInputs_random, 
                                 testNumber=TestStatus.getTestSetup(), 
                                 buildTrees=False, buildFaces=False, buildReference=True,
                                 configFile=configFile)
    os.remove(configFile)
    
def getConfigFile(matchingAlgorithm="greedy"):
    tempConfigFile = getTempFile(rootDir="./", suffix=".xml")
    config = ET.parse(os.path.join(cactusRootPath(), "pipeline", "cactus_workflow_config.xml")).getroot()
    #Set the matching algorithm
    config.find("reference").attrib["matching_algorithm"] = matchingAlgorithm
    #Now print the file..
    fileHandle = open(tempConfigFile, 'w')
    ET.ElementTree(config).write(fileHandle)
    fileHandle.close()
    return tempConfigFile
    
def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
