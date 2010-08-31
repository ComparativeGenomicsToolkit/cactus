import unittest
import sys
import os
import random

from sonLib.bioio import logger
from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import system
from sonLib.bioio import getTempDirectory

from cactus.shared.common import runCactusSetup
from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusWorkflowExperimentConfig
from cactus.shared.test import getCactusDiskDatabaseString

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
            experiment = getCactusWorkflowExperimentConfig(tempDir, sequences, newickTreeString)
            cactusDiskDatabaseString = getCactusDiskDatabaseString(experiment)
            #Make the experiment file
            #experimentFile = makeCactusWorkflowExperimentFile(tempDir, experiment)
            
            runCactusSetup(cactusDiskDatabaseString, sequences, newickTreeString, debug=True)
            runCactusSetup(cactusDiskDatabaseString, sequences, newickTreeString, debug=True)
            system("rm -rf %s" % tempDir)
            logger.info("Finished test %i of cactus_setup.py", test) 
 
def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
