#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys
import os
import xml.etree.ElementTree as ET
import random

from cactus.shared.test import parseCactusSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempFile
from sonLib.bioio import system

from cactus.shared.common import cactusRootPath

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import runWorkflow_multipleExamples

class TestCase(unittest.TestCase):
    def testCactusCore_Random(self):
        for test in xrange(TestStatus.getTestSetup()):
            randomConfigFile=getRandomConfigFile()
            runWorkflow_multipleExamples(getCactusInputs_random, 
                                         buildTrees=False, buildFaces=False, buildReference=False, configFile=randomConfigFile)
            os.remove(randomConfigFile)
        
    def testCactusCore_Blanchette(self):
        runWorkflow_multipleExamples(getCactusInputs_blanchette, 
                                     testRestrictions=(TestStatus.TEST_SHORT,), inverseTestRestrictions=True, 
                                     buildTrees=False, buildFaces=False, buildReference=False)
        
def getRandomConfigFile():
    tempConfigFile = getTempFile(rootDir="./", suffix=".xml")
    config = ET.parse(os.path.join(cactusRootPath(), "pipeline", "cactus_workflow_config.xml")).getroot()
    #Mess with the number of iterations and the parameters for the iteration..
    iterations = config.find("alignment").find("iterations")
    i = iterations.findall("iteration")
    #Remove all the iterations bar one..
    i.reverse()
    for iteration in i:
        if iteration.attrib["type"] == "blast":
            for iteration2 in i:
                if iteration2 != iteration:
                    iterations.remove(iteration2)
            break
    #Now make random parameters..
    iteration.attrib["number"] = "0"
    core = iteration.find("core")
    annealingRounds = 1 + int(random.random() * 10)
    
    core.attrib["annealingRounds"] = " ".join([ str(1 + int(random.random() * 10)) for i in xrange(annealingRounds) ])
    deannealingRounds = list(set([ 1 + int(random.random() * 10) for i in xrange(int(random.random() * 10)) ]))
    deannealingRounds.sort()
    core.attrib["deannealingRounds"] = " ".join([ str(i) for i in deannealingRounds ])
    core.attrib["trim"] = " ".join([ str(1 + int(random.random() * 5)) for i in xrange(annealingRounds) ])
    
    core.attrib["alignRepeatsAtLoop"] = str(random.random() * annealingRounds)
    
    core.attrib["minimumTreeCoverage"] = str(random.random())
    core.attrib["blockTrim"] = str(int(random.random() * 5))
    core.attrib["ignoreAllChainsLessThanMinimumTreeCoverage"] = str(random.choice([0, 1]))
    
    #Now print the file..
    fileHandle = open(tempConfigFile, 'w')
    ET.ElementTree(config).write(fileHandle)
    fileHandle.close()
    system("cat %s" % tempConfigFile)
    return tempConfigFile
    
def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
