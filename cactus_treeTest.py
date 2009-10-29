"""Tests the core of the cactus tree pipeline.
"""

import unittest

import os
import sys
import random

from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempDirectory
from sonLib.bioio import getTempFile
from sonLib.bioio import logger
from sonLib.bioio import system

from cactus.cactus_common import runCactusSetup
from cactus.cactus_common import runCactusAligner
from cactus.cactus_common import runCactusCore
from cactus.cactus_common import runCactusTree
from cactus.cactus_common import getRandomCactusInputs
from cactus.cactus_common import runCactusCheckReconstructionTree

class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.testNo = TestStatus.getTestSetup(1, 5, 20, 100)
        self.tempDir = getTempDirectory(os.getcwd())
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        system("rm -rf %s" % self.tempDir)
    
    def testCactusTree(self):
        for test in xrange(self.testNo):
            sequenceDirs, newickTreeString = getRandomCactusInputs(tempDir=getTempDirectory(self.tempDir))
            runPipe(sequenceDirs, newickTreeString, self.tempDir, useDummy=True, writeDebugFiles=True,
                    randomAtomParameters=True)
            
def runPipe(sequenceDirs, newickTreeString, tempDir, useDummy=False, writeDebugFiles=False, randomAtomParameters=False):
    tempAlignmentFile = getTempFile(rootDir=tempDir)
    tempReconstructionDirectory = os.path.join(getTempDirectory(tempDir), "tempReconstruction")
    
    runCactusSetup(tempReconstructionDirectory, sequenceDirs, 
                   newickTreeString, getTempDirectory(tempDir))
    runCactusAligner(tempReconstructionDirectory, tempAlignmentFile,
                     tempDir=getTempDirectory(tempDir), useDummy=useDummy)
    runCactusCore(tempReconstructionDirectory, tempAlignmentFile, 
                    tempDir=getTempDirectory(tempDir), writeDebugFiles=writeDebugFiles)
    runCactusTree(tempReconstructionDirectory, tempDir=getTempDirectory(tempDir))
    
    system("rm -rf %s %s" % (tempReconstructionDirectory, tempAlignmentFile))
    
    logger.info("Ran the test of the reconstruction program okay")
        
def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()