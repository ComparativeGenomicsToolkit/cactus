import unittest
import sys
import os
import xml.etree.ElementTree as ET
import random

from sonLib.bioio import logger
from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import system
from sonLib.bioio import getTempDirectory
from sonLib.bioio import fastaRead

from cactus.cactus_common import runCactusSetup
from cactus.cactus_common import getRandomCactusInputs

class TestCase(unittest.TestCase):

    def setUp(self):
        self.testNo = TestStatus.getTestSetup()
        self.tempFiles = []
        self.tempDir = getTempDirectory()
        self.tempReconstructionDirectory = os.path.join(self.tempDir, "tempReconstruction")
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        for tempFile in self.tempFiles:
            os.remove(tempFile)
        unittest.TestCase.tearDown(self)
        #system("rm -rf %s" % self.tempDir)
        
    def testCactusSetup(self):
        for test in xrange(self.testNo): 
            ##########################################
            #Make random inputs
            ##########################################
            
            sequenceNumber = random.choice(xrange(100))
            sequenceDirs, newickTreeString = getRandomCactusInputs(tempDir=getTempDirectory(self.tempDir), sequenceNumber=sequenceNumber)
            cactusTempDir=getTempDirectory(self.tempDir)
            
            runCactusSetup(self.tempReconstructionDirectory, sequenceDirs, 
                               newickTreeString, cactusTempDir)
            
            system("rm -rf %s" % self.tempReconstructionDirectory)
            
            return
            
            runCactusSetup(self.tempReconstructionDirectory, sequenceDirs, 
                               newickTreeString, cactusTempDir) #Run it twice to check the job is atomic.
    
            ##########################################
            #Evaluate the output.
            ##########################################
            
            #Cat the file to the screen
            logger.info("The top level reconstruction problem")
            system("cat %s" % os.path.join(self.tempReconstructionDirectory, "reconstructionProblem.xml"))
            system("cat %s" % os.path.join(self.tempReconstructionDirectory, "fastaMap.xml"))
    
            #File
            reconProblem = ET.parse(os.path.join(self.tempReconstructionDirectory, "reconstructionProblem.xml")).getroot()
            fastaMap = ET.parse(os.path.join(self.tempReconstructionDirectory, "fastaMap.xml")).getroot()
    
            #Check the fasta sequences have been properly processed.
            assert len(reconProblem.find("sequences").findall("sequence")) == sequenceNumber
            assert len(fastaMap.find("fasta_map").findall("fasta")) == sequenceNumber
            
            fMap = {}
            seqFiles = set()
            for fasta in fastaMap.find("fasta_map").findall("fasta"):
                i = False
                for sequence in reconProblem.find("sequences").findall("sequence"):
                    if sequence.attrib["contig"] == fasta.attrib["contig"]:
                        assert i == False
                        i = True
                        assert fasta.attrib["event"] == sequence.attrib["event"]
                        fMap[fasta.attrib["header"]] = os.path.join(self.tempReconstructionDirectory, sequence.attrib["sequence_file"])
                        assert sequence.attrib["sequence_file"] not in seqFiles
                        seqFiles.add(sequence.attrib["sequence_file"])
                assert i == True
            
            for sequenceDir in sequenceDirs:
                for sequenceFile in os.listdir(sequenceDir):
                    fileHandle = open(os.path.join(sequenceDir, sequenceFile), 'r')
                    for fastaHeader, sequence in fastaRead(fileHandle):
                        assert fMap.has_key(fastaHeader)
                        fileHandle2 = open(fMap[fastaHeader], 'r')
                        sequence2 = fileHandle2.read()
                        assert sequence == sequence2
                        fileHandle2.close()
                    fileHandle.close()
            
            logger.info("\nChecked the sequences in the reconstruction problem")
    
    
            system("rm -rf %s" % self.tempReconstructionDirectory)
            
            logger.info("Finished test of cactus_setup.py")

def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()