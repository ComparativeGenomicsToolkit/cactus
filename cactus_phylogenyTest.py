"""Tests the tree building of the cactus pipeline.
"""

import unittest
import os
import sys

from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempDirectory
from sonLib.bioio import getTempFile
from sonLib.bioio import logger
from sonLib.bioio import system

from cactus.cactus_common import runCactusSetup
from cactus.cactus_common import runCactusAligner
from cactus.cactus_common import runCactusCore
from cactus.cactus_common import runCactusPhylogeny
from cactus.cactus_common import getRandomCactusInputs
from cactus.cactus_common import runCactusCheck

class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.testNo = TestStatus.getTestSetup(1, 5, 20, 100)
        self.tempDir = getTempDirectory(os.getcwd())
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        system("rm -rf %s" % self.tempDir)
        system("rm -rf pinchGraph1.dot pinchGraph2.dot pinchGraph3.dot pinchGraph4.dot cactusGraph1.dot cactusGraph2.dot cactusGraph3.dot net1.dot net2.dot net3.dot pinchGraph5.dot pinchGraph6.dot")
        system("rm -rf pinchGraph1.pdf pinchGraph2.pdf pinchGraph3.pdf pinchGraph4.pdf cactusGraph1.pdf cactusGraph2.pdf cactusGraph3.pdf net1.pdf net2.pdf net3.pdf pinchGraph5.pdf pinchGraph6.pdf")
    
    def testCactusTree(self):
        for test in xrange(self.testNo):
            sequenceDirs, newickTreeString = getRandomCactusInputs(tempDir=getTempDirectory(self.tempDir))
            runPipe(sequenceDirs, newickTreeString, self.tempDir, useDummy=True, writeDebugFiles=False,
                    randomAtomParameters=True)
            
def runPipe(sequenceDirs, newickTreeString, tempDir, useDummy=False, writeDebugFiles=False, randomAtomParameters=False):
    tempAlignmentFile = getTempFile(rootDir=tempDir)
    tempReconstructionDirectory = os.path.join(getTempDirectory(tempDir), "tempReconstruction")
    
    runCactusSetup(tempReconstructionDirectory, sequenceDirs, 
                   newickTreeString, getTempDirectory(tempDir))
    runCactusAligner(tempReconstructionDirectory, tempAlignmentFile,
                     tempDir=getTempDirectory(tempDir), useDummy=useDummy)
    runCactusCore(tempReconstructionDirectory, tempAlignmentFile, 
                    tempDir=getTempDirectory(tempDir), writeDebugFiles=writeDebugFiles,
                    maximumEdgeDegree=10000,
                    proportionOfAtomsToKeep=1.0,
                    discardRatio=0.0,
                    minimumTreeCoverage=0.1,
                    minimumChainLength=1)
    runCactusPhylogeny(tempReconstructionDirectory, tempDir=getTempDirectory(tempDir))
    runCactusCheck(tempReconstructionDirectory)
    
    system("rm -rf %s %s" % (tempReconstructionDirectory, tempAlignmentFile))
    
    logger.info("Ran the test of the reconstruction program okay")
        
def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()