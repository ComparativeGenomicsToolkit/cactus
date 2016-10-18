#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys
import os
import random

from sonLib.bioio import logger
from sonLib.bioio import TestStatus
from sonLib.bioio import system
from sonLib.bioio import getTempDirectory

from cactus.shared.common import runCactusSetup
from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusWorkflowExperimentForTest

class TestCase(unittest.TestCase):

    def setUp(self):
        self.testNo = TestStatus.getTestSetup()
        unittest.TestCase.setUp(self)
       
    def testCactusSetup(self): 
        """Creates a bunch of random inputs and then passes them to cactus setup.
        """
        for test in xrange(self.testNo): 
            tempDir = getTempDirectory(os.getcwd())
            sequenceNumber = random.choice(xrange(100))
            sequences, newickTreeString = getCactusInputs_random(tempDir=tempDir, sequenceNumber=sequenceNumber)
            
            #Setup the flower disk.
            experiment = getCactusWorkflowExperimentForTest(sequences, newickTreeString, tempDir)
            cactusDiskDatabaseString = experiment.getDiskDatabaseString() 
            
            runCactusSetup(cactusDiskDatabaseString, sequences, newickTreeString)
            runCactusSetup(cactusDiskDatabaseString, sequences, newickTreeString)
            
            experiment.cleanupDb()
            system("rm -rf %s" % tempDir)
            logger.info("Finished test %i of cactus_setup.py", test) 
 
def main():
    unittest.main()
        
if __name__ == '__main__':
    main()
