import unittest
import os
import time
from sonLib.bioio import popenPush, system
from sonLib.bioio import TestStatus
from sonLib.bioio import fastaRead
from sonLib.bioio import getTempDirectory

"""Base case used for testing the preprocessor and lastz repeat masking
"""

class TestCase(unittest.TestCase):
    
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.encodeRegion = "ENm001"
        self.encodePath = os.path.join(TestStatus.getPathToDataSets(), "MAY-2005")
        self.regionPath = os.path.join(self.encodePath, self.encodeRegion)
        self.tempDir = getTempDirectory(os.getcwd())
        self.tempOutputFile = os.path.join(self.tempDir, "results1.txt")
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        system("rm -rf %s" % self.tempDir)
        
    def checkSequenceSetsEqualModuloSoftMasking(self, sequences1, sequences2):
        self.assertEquals(sequences1.keys(), sequences2.keys())
        for seqName in sequences1.keys():
            sequence1 = sequences1[seqName]
            sequence2 = sequences2[seqName]
            self.assertEquals(sequence1.upper(), sequence2.upper())
        
def getSequences(sequenceFile):
    sequences = {}
    fileHandle = open(sequenceFile, "r")
    for header, sequence in fastaRead(fileHandle):
        sequences[header] = sequence
    fileHandle.close()
    return sequences

def getMaskedBases(sequences):
    maskedBases = set()
    for header in sequences.keys():
        sequence = sequences[header]
        for i in xrange(len(sequence)):
            base = sequence[i]
            if base.upper() != base or base == 'N':
                maskedBases.add((header, i, base))
    return maskedBases

def getLowerCaseBases(sequenceFile):
    #Counts lower case bases in fasta sequences
    from sonLib.bioio import fastaRead
    totalMasked = 0
    total = 0
    fileHandle = open(sequenceFile, "r")
    for header, sequence in fastaRead(fileHandle):
        for base in sequence:
            if base != base.upper():
                totalMasked += 1
        total += len(sequence)
    fileHandle.close()
    return total, totalMasked
        
if __name__ == '__main__':
    unittest.main()
