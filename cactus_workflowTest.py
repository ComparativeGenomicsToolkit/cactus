import unittest
import os
import sys

from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import runGraphViz

from cactus.cactus_common import getRandomCactusInputs
from cactus.cactus_common import runCactusWorkflow
from cactus.cactus_common import runCactusTreeViewer
from cactus.cactus_common import runCactusBlockGraphViewer
from cactus.cactus_common import runCactusCheck
from cactus.cactus_common import runCactusTreeStats

from workflow.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete

class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.testNo = TestStatus.getTestSetup(1, 100, 0, 0)
        self.sequenceNumber = TestStatus.getTestSetup(5, 50, 50, 100)
        self.tempFiles = []
        self.tempDir = getTempDirectory(os.getcwd())
        self.tempReconstructionDirectory = os.path.join(self.tempDir, "tempReconstruction")
        if os.system("parasol status") == 0:
            self.batchSystem = "parasol"
        else:
            self.batchSystem = "single_machine"
        #self.batchSystem = "single_machine"
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
                               batchSystem=self.batchSystem, buildTrees=True, buildAdjacencies=True)
            runCactusCheck(self.tempReconstructionDirectory)
            #Run the checker to check the file is okay.
            runJobTreeStatusAndFailIfNotComplete(jobTreeDir)
            #Cleanup this test
            system("rm -rf %s %s" % (self.tempReconstructionDirectory, jobTreeDir))
            logger.info("Finished this round of test")

    def testCactusWorkflow_Blanchette(self): 
        """Runs the workflow on blanchette's simulated (colinear) regions.
        """
        if TestStatus.getTestStatus() in (TestStatus.TEST_LONG,):
            blanchettePath = os.path.join(TestStatus.getPathToDataSets(), "blanchettesSimulation")
            
            newickTreeFile = os.path.join(blanchettePath, "tree.newick")
            for region in xrange(0, 1):
                sequences = [ os.path.join(blanchettePath, 
                                               ("%.2i.job" % region), species) \
                                 for species in ("HUMAN", "CHIMP", "BABOON", "MOUSE", "RAT", "DOG", "CAT", "PIG", "COW") ] #Same order as tree
                
                outputDir = os.path.join(TestStatus.getPathToDataSets(), "cactus", 
                                                            "blanchettesRegionsTest", str(region))
                
                runWorkflow(sequences, newickTreeFile, outputDir, self.tempDir, batchSystem=self.batchSystem, 
                            cactusTreeGraphFile=os.path.join(outputDir, "cactusTree.dot"),
                            cactusTreeGraphPDFFile=os.path.join(outputDir, "cactusTree.pdf"),
                            blockGraphFile=os.path.join(outputDir, "blockGraph.dot"),
                            blockGraphPDFFile=os.path.join(outputDir, "blockGraph.pdf"),
                            cactusTreeStatsFile=os.path.join(outputDir, "cactusTreeStats.xml"),
                            buildTrees=False, buildAdjacencies=False)
        
    def testCactusWorkflow_Encode(self): 
        """Run the workflow on the encode pilot regions.
        """
        if TestStatus.getTestStatus() in (TestStatus.TEST_LONG,):
            encodeDatasetPath = os.path.join(TestStatus.getPathToDataSets(), "MAY-2005")
            encodeResultsPath = os.path.join(TestStatus.getPathToDataSets(), "cactus", "encodeRegionsTest")
            newickTreeFile = os.path.join(encodeDatasetPath, "reducedTree.newick")
            
            for encodeRegion in [ "ENm00" + str(i) for i in xrange(1, 5) ]:
                sequences = [ os.path.join(encodeDatasetPath, encodeRegion, ("%s.%s.fa" % (species, encodeRegion))) for\
                             species in ("human", "chimp", "baboon", "mouse", "rat", "dog", "cow") ]
                outputDir = os.path.join(encodeResultsPath, encodeRegion)
                
                runWorkflow(sequences, newickTreeFile, outputDir, self.tempDir, batchSystem=self.batchSystem,
                            cactusTreeStatsFile=os.path.join(outputDir, "cactusTreeStats.xml"))
    
    def testCactusWorkflow_Chromosomes(self):
        #Tests cactus_core on the alignment of 4 whole chromosome X's, human, chimp, mouse, dog.
        if TestStatus.getTestStatus() in (TestStatus.TEST_VERY_LONG,):
            chrXPath = os.path.join(TestStatus.getPathToDataSets(), "chr_x")
            chrXResultsPath = os.path.join(TestStatus.getPathToDataSets(), "cactus", "chrX")
            #sequences = [ os.path.join(chrXPath, seqFile) for seqFile in ("hg18.fa", "mouse_chrX.fa") ]
            sequences = [ os.path.join(chrXPath, seqFile) for seqFile in ("hg18.fa", "panTro2.fa", "mouse_chrX.fa", "dog_chrX.fa") ]
            outputDir = chrXResultsPath
            newickTreeFile = os.path.join(chrXPath, "newickTree.txt")
            #newickTreeString = '(((h:0.006969, c:0.009727):0.13, m:0.36):0.02, d:0.15);'
            runWorkflow(sequences, newickTreeFile, outputDir, self.tempDir, batchSystem=self.batchSystem,
                    cactusTreeStatsFile=os.path.join(outputDir, "chrXStats.xml"))
    
def runWorkflow(sequences, newickTreeFile, outputDir, tempDir, 
                batchSystem="single_machine",
                cactusTreeGraphFile=None, cactusTreeGraphPDFFile=None, 
                blockGraphFile=None, blockGraphPDFFile=None,
                cactusTreeStatsFile=None, buildTrees=False, buildAdjacencies=False):
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
    
    runCactusWorkflow(netDisk, sequences, newickTreeString, jobTreeDir, 
                      batchSystem=batchSystem, buildTrees=buildTrees, buildAdjacencies=buildAdjacencies)
    logger.info("Ran the the workflow")
    
    runJobTreeStatusAndFailIfNotComplete(jobTreeDir)
    logger.info("Checked the job tree dir")
    
    if buildTrees:
        runCactusCheck(netDisk, checkTrees=True)
    else:
        runCactusCheck(netDisk)
    logger.info("Checked the cactus tree")
    
    if cactusTreeGraphFile != None:
        runCactusTreeViewer(cactusTreeGraphFile, netDisk)
        runGraphViz(cactusTreeGraphFile, cactusTreeGraphPDFFile)
        
    if cactusTreeStatsFile != None:
        runCactusTreeStats(netDisk, cactusTreeStatsFile)
    
    """
    if blockGraphFile != None:
        runCactusBlockGraphViewer(blockGraphFile, reconstructionTree)
        runGraphViz(blockGraphFile, blockGraphPDFFile, "neato")
    """
        
    #Cleanup
    system("rm -rf %s" % jobTreeDir)
    logger.info("Successfully ran job for the problem")

def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
