"""Tests the workflow pipeline using jobTree with a repeated crashes of parasol etc.
"""

import unittest
import os
import sys
import time

from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger
from sonLib.bioio import system

from cactus.cactus_common import runCactusCheckReconstructionTree
from cactus.cactus_common import getRandomCactusInputs
from cactus.cactus_common import runCactusWorkflow

from workflow.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete

from workflow.jobTree.jobTreeParasolCrashTest import ParasolAndMasterKiller
from workflow.jobTree.jobTreeParasolCrashTest import parasolStop

class TestCase(unittest.TestCase):
     
    def setUp(self):
        self.testNo = TestStatus.getTestSetup(1, 1, 5, 100)
        self.alignmentIterations = TestStatus.getTestSetup(1, 1, 1, 2)
        self.sequenceNumber = TestStatus.getTestSetup(5, 50, 50, 100)
        self.tempFiles = []
        self.tempDir = getTempDirectory(os.getcwd())
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        for tempFile in self.tempFiles:
            os.remove(tempFile)
        unittest.TestCase.tearDown(self)
        system("rm -rf %s" % self.tempDir)
       
    def testCactusWorkflowWithCrashes(self):
        """Runs the tests across some simulated regions, crashing the system peridodically to see if the workflow recovers.
        """
        for test in xrange(self.testNo): 
            sequenceDirs, newickTreeString = getRandomCactusInputs(tempDir=getTempDirectory(self.tempDir), sequenceNumber=self.sequenceNumber)
            tempReconstructionDirectory = os.path.join(self.tempDir, "tempReconstruction")
            jobTreeDir = os.path.join(getTempDirectory(self.tempDir), "jobTree")
            parasolAndMasterKiller = ParasolAndMasterKiller()
            parasolAndMasterKiller.start()
            while True:
                while True:
                    try:
                        runCactusWorkflow(tempReconstructionDirectory, sequenceDirs, newickTreeString, 
                                          jobTreeDir, logLevel='INFO', batchSystem="parasol", rescueJobFrequency=5, alignmentIterations=self.alignmentIterations)
                    except RuntimeError:
                        logger.info("The job tree master ended with an error, restarting")
                        time.sleep(5)
                        continue
                    logger.info("The job tree master ended, with an okay exit value (using parasol)")
                    break
                try:
                    runJobTreeStatusAndFailIfNotComplete(jobTreeDir) #Check the state of the job files
                except RuntimeError:
                    logger.info("Not completed all the jobs yet")
                    continue
                break
            parasolStop()
            #Run the checker to check the file is okay.
            runJobTreeStatusAndFailIfNotComplete(jobTreeDir)
            runCactusCheckReconstructionTree(tempReconstructionDirectory)
            #Cleanup this test
            system("rm -rf %s %s" % (tempReconstructionDirectory, jobTreeDir))
            logger.info("Finished this round of test")

def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()