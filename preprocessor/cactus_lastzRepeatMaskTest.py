import unittest
import os
from sonLib.bioio import popenPush, system
from sonLib.bioio import TestStatus
from sonLib.bioio import fastaRead
from sonLib.bioio import getTempDirectory

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
        
    def testLastzRepeatMask(self):
        #Demo sequences
        sequenceFiles = [ os.path.join(self.encodePath, self.encodeRegion, "%s.ENm001.fa" % species) for species in 'human', "hedgehog" ]
        #Max occurrences
        maxOccurrences = [ 1, 2, 5, 10, 20 ]
        
        for sequenceFile in sequenceFiles:
            #Parse sequences into dictionary
            originalSequences = getSequences(sequenceFile)
            #Get the masked bases
            maskedBasesOriginal = getMaskedBases(originalSequences)
            #Total bases
            totalBases = sum([ len(i) for i in originalSequences.values() ])
            #Calculate number of hard masked bases
            totalNBases = len([ (header, i, base) for (header, i, base) in maskedBasesOriginal if base.upper() == "N" ])
            
            for maxOccurrence in maxOccurrences:
                
                #Run lastz repeat masker
                command = "cactus_lastzRepeatMask.py --proportionSampled=1.0 --minPeriod=%s --unmask --lastzOpts='--step=1 --ambiguous=iupac --nogapped' %s %s" % \
                       (maxOccurrence, sequenceFile, self.tempOutputFile)
                print "Running command", command
                popenPush(command, sequenceFile)
            
                #Parse lastz masked sequences into dictionary
                lastzSequences = getSequences(self.tempOutputFile)
            
                #Check the sequences are the same modulo masking
                self.assertEquals(originalSequences.keys(), lastzSequences.keys())
                for seqName in originalSequences.keys():
                    originalSequence = originalSequences[seqName]
                    lastzSequence = lastzSequences[seqName]
                    self.assertEquals(originalSequence.upper(), lastzSequence.upper())
            
                #Compare the proportion of bases masked by lastz with original repeat masking
                maskedBasesOriginal = getMaskedBases(originalSequences)
                maskedBasesLastzMasked = getMaskedBases(lastzSequences)
                
                print " For the sequence file ", sequenceFile, \
                 " the total number of sequences is ", len(originalSequences), \
                 " the total number of bases ", totalBases, \
                 " the number of bases originally masked was: ", len(maskedBasesOriginal),\
                 " the number of bases masked after running lastz repeat masking is: ", len(maskedBasesLastzMasked), \
                 " the intersection of these masked sets is: ", len(maskedBasesLastzMasked.intersection(maskedBasesOriginal)), \
                 " the total number of bases that are Ns ", totalNBases, \
                 " lastz was filter for max-occurrences of more than : ", maxOccurrence

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
        
if __name__ == '__main__':
    unittest.main()