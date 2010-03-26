import unittest
import os
import sys
import random

from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import fastaAlignmentRead
from sonLib.bioio import fastaWrite
from sonLib.bioio import fastaAlignmentWrite
from sonLib.bioio import fastaReadHeaders

from cactus.shared.cactus_common import runCactusWorkflow

from workflow.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete

"""Using this set of tests I have found that:

Spanning tree no / Cut off / Sens / Spec

1    0.9    0.5     0.95
2    0.9    0.666394687629 0.959811182484
2    0.85   0.695476098788 0.936469685867  
2    0.8    0.714650750468 0.901942148024 
2    0.75   0.731373533136 0.857005520637
2    0.7    0.747936598984 0.820758253164
2    0.65    0.76194419572 0.77428747256
3    0.9    0.71    0.89
10    0.9    0.84    0.86

2 0.5 MaxEdge 0.638431294206 0.959093882956
2 0.6 MaxEdge 0.667422055278 0.966186684277
2 0.65 MaxEdge 0.671937292932 0.971078624611
2 0.7 MaxEdge 0.679474200144 0.975266205164
2 0.75 MaxEdge 0.676702897492 0.978594642833
2 0.8 MaxEdge 0.675537396377 0.98162754782

All pairs 0.9 0.87 0.81

By default we us 2/0.85

"""

class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.testNo = TestStatus.getTestSetup(5, 100, 0, 0)
        self.tempFiles = []
        self.tempDir = getTempDirectory(os.getcwd())
        self.tempReconstructionDirectory = os.path.join(self.tempDir, "tempReconstruction")
        self.batchSystem = "parasol"
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        for tempFile in self.tempFiles:
            os.remove(tempFile)
        unittest.TestCase.tearDown(self)
        system("rm -rf %s" % self.tempDir)

    def testCactusWorkflow_Blanchette(self): 
        """Runs the workflow on blanchette's simulated (colinear) regions.
        """
        trueAlignment = os.path.join(TestStatus.getPathToDataSets(), "blanchettesSimulation", "00.job", "true.mfa")
        
        #Load the true alignment.
        columnAlignment = [ i for i in  fastaAlignmentRead(trueAlignment) ]
        fastaHeaders = [ i for i in fastaReadHeaders(trueAlignment) ]
        sequenceNumber = 9
        
        #The tree
        newickTreeString = "((((HUMAN:0.006969, CHIMP:0.009727):0.025291, BABOON:0.044568):0.11,(RAT:0.072818, MOUSE:0.081244):0.260342):0.023260,((DOG:0.07, CAT:0.07):0.087381,(PIG:0.06, COW:0.06):0.104728):0.04);"
        
        parentResultsFile = os.path.join(self.tempDir, "finalResults.xml")
        
        random.seed(1)
        
        for test in xrange(self.testNo):
            logger.info("Starting the test %i" % test)
            #Get random dir
            testDir = getTempDirectory(self.tempDir)
            
            #Choose random sub-range of alignment.
            columnStart = random.choice(xrange(len(columnAlignment)))
            subAlignment = columnAlignment[columnStart:columnStart+random.choice(xrange(500))]
            logger.info("Got a sub alignment, it is %i columns long" % len(subAlignment))
            
            #Output the 'TRUE' alignment file
            trueMFAFile = os.path.join(testDir, "true.mfa")
            fastaAlignmentWrite(subAlignment, fastaHeaders, len(fastaHeaders), trueMFAFile)
            trueMAFFile = os.path.join(testDir, "true.maf")
            system("eval_MFAToMAF --mFAFile %s --outputFile %s" % (trueMFAFile, trueMAFFile))
            system("cat %s" % trueMAFFile)
            
            #Get sequences
            sequences = [ (fastaHeaders[seqNo], "".join([ column[seqNo] for column in subAlignment if column[seqNo] != '-' ])) for seqNo in xrange(sequenceNumber) ]
            logger.info("Got the sequences")
            
            #Write sequences into temp files
            tempFastaFiles = []
            for seqNo in xrange(sequenceNumber):
                header, sequence = sequences[seqNo]
                logger.info("Making temp file for header: %s, seq: %s" % (header, sequence))
                tempFastaFile = os.path.join(testDir, "%i.fa" % seqNo)
                tempFastaFiles.append(tempFastaFile)
                fileHandle = open(tempFastaFile, "w")
                fastaWrite(fileHandle, header, sequence)
                fileHandle.close()
            logger.info("Got the temp sequence files")
            
            #Run the workflow
            netDisk = os.path.join(testDir, "netDisk")
            jobTree = os.path.join(testDir, "jobTree")
            runCactusWorkflow(netDisk, tempFastaFiles, newickTreeString, jobTree)
            logger.info("Ran the the workflow")
            
            #Check the output alignment
            runJobTreeStatusAndFailIfNotComplete(jobTree)
            logger.info("Checked the job tree dir")
            
            #Now get mafs for the region.
            mAFFile = os.path.join(testDir, "net.maf")
            system("cactus_MAFGenerator --netName 0 --netDisk %s --outputFile %s" % (netDisk, mAFFile))
            logger.info("Got the MAFs from the net disk")
            system("cat %s" % mAFFile)
            
            statsFile = os.path.join(testDir, "stats.xml")
            system("cactus_treeStats --netDisk %s --netName 0 --outputFile %s" % (netDisk, statsFile))
            system("cat %s" % statsFile)
            logger.info("Got the cactus tree stats")
            
            #Now compare the mafs to the output.
            resultsFile = os.path.join(testDir, "results.xml")
            system("eval_MAFComparator --mAFFile1 %s --mAFFile2 %s --outputFile %s" % (trueMAFFile, mAFFile, resultsFile))
            logger.info("Ran the maf comparator")
            
            system("cat %s" % resultsFile)
            
            if test > 0:
                system("eval_mergeMAFComparatorResults.py --results1 %s --results2 %s --outputFile %s" % (parentResultsFile, resultsFile, parentResultsFile))
                logger.info("Merging the results")
            else:
                system("cp %s %s" % (resultsFile, parentResultsFile))
                logger.info("Copying the results")
            
            #Cleanup
            system("rm -rf %s" % testDir)
            logger.info("Successfully ran test for the problem")
        
        system("cat %s" % parentResultsFile)
        logger.info("The final results file")

def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
