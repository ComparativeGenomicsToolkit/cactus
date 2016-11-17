#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import os
import sys
import random

from sonLib.bioio import TestStatus
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import fastaAlignmentRead
from sonLib.bioio import fastaWrite
from sonLib.bioio import fastaAlignmentWrite
from sonLib.bioio import fastaReadHeaders
from sonLib.bioio import getLogLevelString

from cactus.shared.common import runCactusWorkflow

from cactus.shared.test import getCactusWorkflowExperimentForTest
from cactus.shared.test import silentOnSuccess
from cactus.shared.common import runToilStatusAndFailIfNotComplete

"""Tests cactus_bar. Requires the installation of cactusTools and mafTools.
"""

class TestCase(unittest.TestCase):
    
    def setUp(self):
        self.testNo = TestStatus.getTestSetup(3, 10, 0, 0)
        self.batchSystem = "parasol"
        unittest.TestCase.setUp(self)

    @unittest.skip("")
    def testCactusWorkflow_Blanchette(self): 
        """Runs the workflow on blanchette's simulated (colinear) regions.
        """
        if "SON_TRACE_DATASETS" not in os.environ:
            return
        for test in xrange(self.testNo):
            tempFiles = []
            tempDir = getTempDirectory(os.getcwd())
            
            trueAlignment = os.path.join(TestStatus.getPathToDataSets(), "blanchettesSimulation", "00.job", "true.mfa")
            
            #Load the true alignment.
            columnAlignment = [ i for i in  fastaAlignmentRead(trueAlignment) ]
            fastaHeaders = [ i for i in fastaReadHeaders(trueAlignment) ]
            sequenceNumber = 9
            
            #The tree
            newickTreeString = "((((HUMAN:0.006969, CHIMP:0.009727):0.025291, BABOON:0.044568):0.11,(RAT:0.072818, MOUSE:0.081244):0.260342):0.023260,((DOG:0.07, CAT:0.07):0.087381,(PIG:0.06, COW:0.06):0.104728):0.04);"
            
            #Get random dir
            testDir = getTempDirectory(tempDir)
            
            #random alignment
            alignmentLength = 5000
            randomStart = random.choice(xrange(len(columnAlignment)-alignmentLength))
            subAlignment = columnAlignment[randomStart:randomStart+alignmentLength]
            logger.info("Got a sub alignment, it is %i columns long" % len(subAlignment))
            
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

            experiment = getCactusWorkflowExperimentForTest(tempFastaFiles, newickTreeString, testDir)
            experimentFile = os.path.join(testDir, "experiment.xml")
            experiment.writeXML(experimentFile)
            cactusDiskDatabaseString = experiment.getDiskDatabaseString()

            toilDir = os.path.join(testDir, "jobStore")

            runCactusWorkflow(experimentFile, toilDir)
            logger.info("Ran the the workflow")

            runToilStatusAndFailIfNotComplete(toilDir)
            logger.info("Checked the toil status")

            #Output the 'TRUE' alignment file
            if os.system("mfaToMaf --help > /dev/null 2>&1") == 0 and\
               os.system("cactus_MAFGenerator --help > /dev/null 2>&1") == 0 and\
               os.system("mafComparator --help > /dev/null 2>&1") == 0 and\
               os.system("cactus_treeStats --help > /dev/null 2>&1") == 0:
                trueMFAFile = os.path.join(testDir, "true.mfa")
                fastaAlignmentWrite(subAlignment, fastaHeaders, len(fastaHeaders), trueMFAFile)
                trueMAFFile = os.path.join(testDir, "true.maf")
                system("mfaToMaf --mfaFile %s --outputFile %s --logLevel %s" % (trueMFAFile, trueMAFFile, getLogLevelString()))
                system("cat %s" % trueMAFFile)
                
                #Now get mafs for the region.
                mAFFile = os.path.join(testDir, "flower.maf")
                system("cactus_MAFGenerator --flowerName 0 --cactusDisk '%s' --outputFile %s --logLevel %s" % (cactusDiskDatabaseString, mAFFile, getLogLevelString()))
                logger.info("Got the MAFs from the flower disk")
                system("cat %s" % mAFFile)
                
                statsFile = os.path.join(testDir, "stats.xml")
                system("cactus_treeStats --cactusDisk '%s' --flowerName 0 --outputFile %s --logLevel %s" % (cactusDiskDatabaseString, statsFile, getLogLevelString()))
                system("cat %s" % statsFile)
                logger.info("Got the cactus tree stats")
                
                #Now compare the mafs to the output.
                resultsFile = os.path.join(testDir, "results.xml")
                system("mafComparator --mafFile1 %s --mafFile2 %s --outputFile %s --logLevel %s" % (trueMAFFile, mAFFile, resultsFile, getLogLevelString()))
                logger.info("Ran the maf comparator")
                
                system("cat %s" % resultsFile)
                
                #Cleanup
                experiment.cleanupDb()
                system("rm -rf %s" % testDir)
                logger.info("Successfully ran test for the problem")
                
            for tempFile in tempFiles:
                os.remove(tempFile)
            system("rm -rf %s" % tempDir)
    
    @silentOnSuccess
    def testPosetAlignerAPI(self):
        """Run all the cactus base aligner CuTests, fail if any of them fail.
        """
        system("cactus_barTests %s" % getLogLevelString())

if __name__ == '__main__':
    unittest.main()
