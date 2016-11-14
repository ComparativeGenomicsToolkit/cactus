#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys
import os
import xml.etree.ElementTree as ET
import random

from sonLib.bioio import TestStatus
from sonLib.bioio import getTempFile
from sonLib.bioio import system
from sonLib.bioio import getLogLevelString

from cactus.shared.common import cactusRootPath

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import runWorkflow_multipleExamples

class TestCase(unittest.TestCase):
    def testCactusCore_Random(self):
        for test in xrange(TestStatus.getTestSetup()):
            randomConfigFile=getRandomConfigFile()
            runWorkflow_multipleExamples(getCactusInputs_random, 
                                         configFile=randomConfigFile)
            os.remove(randomConfigFile)
        
    def testCactusCore_Blanchette(self):
        runWorkflow_multipleExamples(getCactusInputs_blanchette, 
                                     testRestrictions=(TestStatus.TEST_SHORT,), inverseTestRestrictions=True)
        
    def testCuTest(self):
        system("stCafTests %s" % getLogLevelString())       
        
def getRandomConfigFile():
    tempConfigFile = getTempFile(rootDir="./", suffix=".xml")
    config = ET.parse(os.path.join(cactusRootPath(), "cactus_config.xml")).getroot()
    cafNode = config.find("caf")
    assert len(config.findall("caf")) == 1
    
    annealingRounds = 1 + int(random.random() * 10)
    cafNode.attrib["annealingRounds"] = " ".join([ str(1 + int(random.random() * 10)) for i in xrange(annealingRounds) ])
    deannealingRounds = list(set([ 1 + int(random.random() * 10) for i in xrange(int(random.random() * 10)) ]))
    deannealingRounds.sort()
    cafNode.attrib["deannealingRounds"] = " ".join([ str(i) for i in deannealingRounds ])
    cafNode.attrib["trim"] = " ".join([ str(1 + int(random.random() * 5)) for i in xrange(annealingRounds) ])
    
    cafNode.attrib["alignRepeatsAtLoop"] = str(random.random() * annealingRounds)
    
    cafNode.attrib["minimumTreeCoverage"] = str(random.random())
    cafNode.attrib["blockTrim"] = str(int(random.random() * 5))
    cafNode.attrib["ignoreAllChainsLessThanMinimumTreeCoverage"] = str(random.choice([0, 1]))
    cafNode.attrib["minimumBlockDegree"] = str(random.choice([0, 5]))
    
    checkNode = config.find("check")
    checkNode.attrib["runCheck"] = "1"
    
    checkNode = config.find("normal")
    checkNode.attrib["iterations"] = "2"
    
    #Now print the file..
    fileHandle = open(tempConfigFile, 'w')
    ET.ElementTree(config).write(fileHandle)
    fileHandle.close()
    if getLogLevelString() == "DEBUG":
        system("cat %s" % tempConfigFile)
    return tempConfigFile
    
if __name__ == '__main__':
    unittest.main()
