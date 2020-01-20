import os
import pytest
import shutil
import unittest
from sonLib.bioio import TestStatus
from sonLib.bioio import fastaRead
from sonLib.bioio import getTempDirectory

from toil.job import Job

"""Base case used for testing the preprocessor and lastz repeat masking
"""

@pytest.mark.blast
@TestStatus.needsTestData
class TestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.encodeRegion = "ENm001"
        self.encodePath = os.path.join(TestStatus.getPathToDataSets(), "MAY-2005")
        self.regionPath = os.path.join(self.encodePath, self.encodeRegion)
        self.tempDir = getTempDirectory(os.getcwd())
        self.tempOutputFile = os.path.join(self.tempDir, "results1.txt")
        self.toilDir = os.path.join(self.tempDir, "toil")
        self.toilOptions = Job.Runner.getDefaultOptions(self.toilDir)
        self.toilOptions.disableCaching = True

    def tearDown(self):
        unittest.TestCase.tearDown(self)
        shutil.rmtree(self.tempDir)

    def checkSequenceSetsEqualModuloSoftMasking(self, sequences1, sequences2):
        self.assertEqual(list(sequences1.keys()), list(sequences2.keys()))
        for seqName in list(sequences1.keys()):
            sequence1 = sequences1[seqName]
            sequence2 = sequences2[seqName]
            self.assertEqual(sequence1.upper(), sequence2.upper())

def getSequences(sequenceFile):
    sequences = {}
    fileHandle = open(sequenceFile, "r")
    for header, sequence in fastaRead(fileHandle):
        sequences[header] = sequence
    fileHandle.close()
    return sequences

def getMaskedBases(sequences):
    maskedBases = set()
    for header in list(sequences.keys()):
        sequence = sequences[header]
        for i in range(len(sequence)):
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
