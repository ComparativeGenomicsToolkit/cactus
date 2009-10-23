import unittest
import os
import sys

from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import runGraphViz

from cactus.cactus_common import runCactusCheckReconstructionTree
from cactus.cactus_common import getRandomCactusInputs
from cactus.cactus_common import runCactusWorkflow
from cactus.cactus_common import runCactusReconstructionTreeViewer
from cactus.cactus_common import runCactusAtomGraphViewer

from workflow.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete
from workflow.jobTree.jobTreeTest import plotJobTreeGraph

class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.testNo = TestStatus.getTestSetup(1, 1, 5, 1)
        self.alignmentIterations = TestStatus.getTestSetup(2, 2, 2, 3)
        self.sequenceNumber = TestStatus.getTestSetup(5, 50, 50, 100)
        self.tempFiles = []
        self.tempDir = getTempDirectory(os.getcwd())
        self.tempReconstructionDirectory = os.path.join(self.tempDir, "tempReconstruction")
        if os.system("parasol status") == 0:
            self.batchSystem = "parasol"
        else:
            self.batchSystem = "single_machine"
        self.batchSystem = "single_machine"
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        for tempFile in self.tempFiles:
            os.remove(tempFile)
        unittest.TestCase.tearDown(self)
        system("rm -rf %s" % self.tempDir)
        
    def testCactusWorkflow_Random(self):
        """Runs the tests across some simulated regions.
        """
        for test in xrange(self.testNo): 
            sequenceDirs, newickTreeString = getRandomCactusInputs(tempDir=getTempDirectory(self.tempDir), sequenceNumber=self.sequenceNumber)
            jobTreeDir = os.path.join(getTempDirectory(self.tempDir), "jobTree")
            runCactusWorkflow(self.tempReconstructionDirectory, sequenceDirs, newickTreeString, jobTreeDir, logLevel='INFO',
                               batchSystem=self.batchSystem, alignmentIterations=self.alignmentIterations)
            #Run the checker to check the file is okay.
            runJobTreeStatusAndFailIfNotComplete(jobTreeDir)
            #runCactusCheckReconstructionTree(self.tempReconstructionDirectory)
            #Cleanup this test
            system("rm -rf %s %s" % (self.tempReconstructionDirectory, jobTreeDir))
            logger.info("Finished this round of test")

    def testCactusWorkflow_Blanchette(self): 
        """Runs the workflow on blanchette's simulated (colinear) regions.
        """
        if TestStatus.getTestStatus() in (TestStatus.TEST_LONG, TestStatus.TEST_VERY_LONG):
            blanchettePath = os.path.join(TestStatus.getPathToDataSets(), "blanchettesSimulation")
            
            newickTreeFile = os.path.join(blanchettePath, "tree.newick")
            for region in xrange(0, 1):
                sequences = [ os.path.join(blanchettePath, 
                                               ("%.2i.job" % region), species) \
                                 for species in ("HUMAN", "CHIMP", "BABOON", "MOUSE", "RAT", "DOG", "CAT", "PIG", "COW") ] #Same order as tree
                
                outputDir = os.path.join(TestStatus.getPathToDataSets(), "cactus", 
                                                            "blanchettesRegionsTest", str(region))
                
                runWorkflow(sequences, newickTreeFile, outputDir, self.tempDir, batchSystem=self.batchSystem, 
                            reconstructionTreeGraphFile=os.path.join(outputDir, "reconTree.dot"),
                            reconstructionTreeGraphPDFFile=os.path.join(outputDir, "reconTree.pdf"),
                            atomGraphFile=os.path.join(outputDir, "atomGraph.dot"),
                            atomGraphPDFFile=os.path.join(outputDir, "atomGraph.pdf"),
                            jobTreeGraphPDFFile=os.path.join(outputDir, "jobTreeGraph.pdf"))
        
    def testCactusWorkflow_Encode(self): 
        """Run the workflow on the encode pilot regions.
        """
        if TestStatus.getTestStatus() in (TestStatus.TEST_LONG, TestStatus.TEST_VERY_LONG):
            encodeDatasetPath = os.path.join(TestStatus.getPathToDataSets(), "MAY-2005")
            encodeResultsPath = os.path.join(TestStatus.getPathToDataSets(), "cactus", "encodeRegionsTest")
            newickTreeFile = os.path.join(encodeDatasetPath, "reducedTree.newick")
            
            for encodeRegion in [ "ENm00" + str(i) for i in xrange(1, 3) ]:
                sequences = [ os.path.join(encodeDatasetPath, encodeRegion, ("%s.%s.fa" % species, encodeRegion)) for\
                             species in ("human", "chimp", "baboon", "mouse", "rat", "dog", "cow") ]
                outputDir = os.path.join(encodeResultsPath, encodeRegion)
                
                runWorkflow(sequences, newickTreeFile, outputDir, self.tempDir, batchSystem=self.batchSystem)

def runWorkflow(sequences, newickTreeFile, outputDir, tempDir, 
                batchSystem="single_machine",
                reconstructionTreeGraphFile=None, reconstructionTreeGraphPDFFile=None, 
                atomGraphFile=None, atomGraphPDFFile=None, 
                jobTreeGraphPDFFile=None):
    fileHandle = open(newickTreeFile, 'r')
    newickTreeString = fileHandle.readline()
    fileHandle.close()
    logger.info("Passed the tree file: %s" % newickTreeString)
    
    if not os.path.isdir(outputDir):
        os.mkdir(outputDir)
    netDisk = os.path.join(outputDir, "netDisk")
    system("rm -rf %s" % netDisk)
    logger.info("Cleaned up any previous net disk: %s" % netDisk)
    
    jobTreeDir = os.path.join(getTempDirectory(tempDir), "jobTree")
    logger.info("Got a job tree dir for the test: %s" % jobTreeDir)
    
    runCactusWorkflow(netDisk, sequences, newickTreeString, jobTreeDir, batchSystem=batchSystem, alignmentIterations=1)
    logger.info("Ran the the workflow")
    
    runJobTreeStatusAndFailIfNotComplete(jobTreeDir)
    logger.info("Checked the job tree dir")
    
    #runCactusCheckReconstructionTree(reconstructionTree)
    #logger.info("Checked the reconstruction tree")
    
    if reconstructionTreeGraphFile != None:
        runCactusReconstructionTreeViewer(reconstructionTreeGraphFile, reconstructionTree)
        runGraphViz(reconstructionTreeGraphFile, reconstructionTreeGraphPDFFile)
     
    if atomGraphFile != None:
        runCactusAtomGraphViewer(atomGraphFile, reconstructionTree)
        runGraphViz(atomGraphFile, atomGraphPDFFile, "neato")
     
    if jobTreeGraphPDFFile != None:
        plotJobTreeGraph(jobTreeDir, jobTreeGraphPDFFile)
        
    #Cleanup
    system("rm -rf %s" % jobTreeDir)
    logger.info("Successfully ran job for the blanchette sequences")

def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()