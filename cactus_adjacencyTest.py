"""Tests the core of the cactus atomiser/reconstruction pipeline.
"""

import unittest

import os
import sys
import xml.etree.ElementTree as ET

from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempDirectory
from sonLib.bioio import getTempFile
from sonLib.bioio import logger
from sonLib.bioio import system

from cactus.cactus_common import runCactusSetup
from cactus.cactus_common import runCactusAligner
from cactus.cactus_common import runCactusCore
from cactus.cactus_common import runCactusAdjacencyBuilder
from cactus.cactus_common import runCactusCheckReconstructionTree
from cactus.cactus_common import getRandomCactusInputs

class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.testNo = TestStatus.getTestSetup(1, 5, 10, 100)
        self.sequenceNumber = TestStatus.getTestSetup(5, 50, 50, 100)
        self.tempFiles = []
        self.tempDir = getTempDirectory()
        self.tempReconstructionDirectory = self.tempDir + "/tempReconstruction"
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        for tempFile in self.tempFiles:
            os.remove(tempFile)
        unittest.TestCase.tearDown(self)
        system("rm -rf %s" % self.tempDir)
        
    def testAdjacencyBuilding(self):
        tempAlignmentFile = getTempFile()
        self.tempFiles.append(tempAlignmentFile)
        
        for test in xrange(self.testNo): 
            ##########################################
            #Make inputs
            ##########################################
            
            sequenceDirs, newickTreeString = getRandomCactusInputs(tempDir=getTempDirectory(self.tempDir), sequenceNumber=self.sequenceNumber)
            runCactusSetup(self.tempReconstructionDirectory, sequenceDirs, 
                           newickTreeString, tempDir=getTempDirectory(self.tempDir))
            runCactusAligner(self.tempReconstructionDirectory, tempAlignmentFile, 
                             tempDir=getTempDirectory(self.tempDir))
            runCactusCore(self.tempReconstructionDirectory, tempAlignmentFile, 
                          tempDir=getTempDirectory(self.tempDir))
        
            logger.info("Setup the reconstruction problems")
            
            def fn(reconstructionProblem):
                tree = ET.parse(os.path.join(self.tempReconstructionDirectory,
                                             reconstructionProblem)).getroot()
                adjacencyComponent = tree.find("adjacency_components")
                adjacencyComponents = adjacencyComponent.findall("adjacency_component")
                for adjacencyComponent in adjacencyComponents:
                    fn(adjacencyComponent.attrib["child_file"])
                runCactusAdjacencyBuilder(self.tempReconstructionDirectory, 
                                          reconstructionProblem, tempDir=getTempDirectory(self.tempDir))
            fn("reconstructionProblem.xml")
            
            logger.info("Ran the cactus adjacency builder okay")
            
            runCactusCheckReconstructionTree(self.tempReconstructionDirectory)
            
            #Run the checker to check the file is okay.
            system("rm -rf %s" % self.tempReconstructionDirectory)
        
def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()