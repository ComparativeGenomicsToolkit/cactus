"""Tests the core of the cactus blockiser/reconstruction pipeline.
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

from cactus.shared.common import runCactusSetup
from cactus.shared.common import runCactusAligner
from cactus.shared.common import runCactusCore
from cactus.shared.common import runCactusCheck
from cactus.shared.common import runCactusExtendNets

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import getCactusInputs_chromosomeX
from cactus.shared.test import getCactusInputs_encode

class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.testNo = TestStatus.getTestSetup(1, 100, 0, 0)
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        system("rm -rf pinchGraph* cactusGraph*")
        
    def testCactusCore_Chromosomes(self):
        """Tests cactus_core on the alignment of 4 whole chromosome X's, human, chimp, mouse, dog.
        """
        if TestStatus.getTestStatus() in (TestStatus.TEST_VERY_LONG,):
            tempDir = getTempDirectory(os.getcwd())
            sequences, newickTreeString = getCactusInputs_chromosomeX()
            runPipe(sequences, newickTreeString, tempDir)
            system("rm -rf %s" % tempDir)
    
    def testCactusCore_Encode(self): 
        """Run the cactus_core on the CFTR encode pilot region.
        """
        if TestStatus.getTestStatus() in (TestStatus.TEST_LONG,):
            tempDir = getTempDirectory(os.getcwd())
            sequences, newickTreeString = getCactusInputs_encode()
            runPipe(sequences, newickTreeString, tempDir)
            system("rm -rf %s" % tempDir)
    
    def testCactusCore_Blanchette(self): 
        """Runs cactus_core on blanchette's simulated (colinear) regions.
        """
        if TestStatus.getTestStatus() in (TestStatus.TEST_MEDIUM,):
            tempDir = getTempDirectory(os.getcwd())
            sequences, newickTreeString = getCactusInputs_blanchette()
            runPipe(sequences, newickTreeString, tempDir)
            system("rm -rf %s" % tempDir)
            
    def testCactusCore_Random(self):
        """Tests core on some random regions.
        """
        for test in xrange(self.testNo):
            tempDir = getTempDirectory(os.getcwd())
            sequences, newickTreeString = getCactusInputs_random(tempDir)
            runPipe(sequences, newickTreeString, tempDir, 
                    useDummy=True, writeDebugFiles=True,
                    randomBlockParameters=True)
            
            ##########################################
            #Make neatos
            ##########################################
            """
            system("neato pinchGraph1.dot -Tpdf > pinchGraph1.pdf") 
            system("neato pinchGraph2.dot -Tpdf > pinchGraph2.pdf")
            system("neato pinchGraph3.dot -Tpdf > pinchGraph3.pdf")
            system("neato cactusGraph.dot -Tpdf > cactusGraph1.pdf")
            """
            system("rm -rf %s" % tempDir)
            
def runPipe(sequenceDirs, newickTreeString, tempDir, useDummy=False, writeDebugFiles=False, randomBlockParameters=False):
    #Redo with flexible parameters...
    tempAlignmentFile = getTempFile(rootDir=tempDir)
    netDisk = os.path.join(getTempDirectory(tempDir), "tempReconstruction")
    
    runCactusSetup(netDisk, sequenceDirs, 
                   newickTreeString)
    
    l = runCactusExtendNets(netDisk, 0, getTempDirectory(tempDir))
                                                  
    if len(l) > 0:
        childNetName, childNetSize = l[0]
    
        runCactusAligner(netDisk, tempAlignmentFile,
                         netName=childNetName,
                         tempDir=getTempDirectory(tempDir), useDummy=useDummy)
        
        system("cat %s" % tempAlignmentFile)
        
        
        logger.info("Constructed the alignments")
        if randomBlockParameters:
            alignUndoLoops = 1 + int(random.random() * 10)
            runCactusCore(netDisk, tempAlignmentFile, 
                          netName=childNetName,
                          writeDebugFiles=writeDebugFiles,
                          alignUndoLoops=alignUndoLoops,
                          maximumEdgeDegree=1+random.random()*10,
                          minimumTreeCoverage=random.random(),
                          minimumBlockLength=random.random()*5,
                          minimumChainLength=random.random()*10,
                          trim=random.random()*5,
                          alignRepeatsAtLoop=int(random.random() * alignUndoLoops))
        else:
            runCactusCore(netDisk, tempAlignmentFile, netName=childNetName,
                          writeDebugFiles=writeDebugFiles)
  
    runCactusCheck(netDisk, recursive=True)
  
    system("rm -rf %s %s" % (netDisk, tempAlignmentFile))
    
    logger.info("Ran the test of the reconstruction program okay")
        
def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
