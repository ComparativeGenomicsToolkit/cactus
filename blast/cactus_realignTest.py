import unittest

import os, sys
from sonLib.bioio import cigarReadFromString, cigarWrite
from sonLib.bioio import popenCatch
from sonLib.bioio import logger
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempFile
from sonLib.bioio import system

class TestCase(unittest.TestCase):
    
    def setUp(self):
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        
    def testCactusRealignDummy(self):
        """Runs cactus realign using the "rescoreOriginalAlignment" mode
        and checks the output is equivalent to what you'd get by just running lastz.
        """
        for seqFile1, seqFile2 in seqFilePairGenerator():
            realignCommand, lastzCommand = getCommands(seqFile1, seqFile2, "--rescoreOriginalAlignment")
            for realignLine, lastzLine in zip([ i for i in popenCatch(realignCommand).split("\n") if i != '' ], 
                                              [ i for i in popenCatch(lastzCommand).split("\n") if i != '' ]):
                realignCigar = cigarReadFromString(realignLine)
                lastzCigar = cigarReadFromString(lastzLine)
                self.assertTrue(realignCigar != None)
                self.assertTrue(realignCigar == lastzCigar)
    
    def testCactusRealign(self):
        """Runs cactus realign using the default parameters and checks that the realigned output cigars align
        the same subsequences.
        """
        for seqFile1, seqFile2 in seqFilePairGenerator():
            realignCommand, lastzCommand = getCommands(seqFile1, seqFile2)
            for realignLine, lastzLine in zip([ i for i in popenCatch(realignCommand).split("\n") if i != '' ], 
                                              [ i for i in popenCatch(lastzCommand).split("\n") if i != '' ]):
                realignCigar = cigarReadFromString(realignLine)
                lastzCigar = cigarReadFromString(lastzLine)
                self.assertTrue(realignCigar.sameCoordinates(lastzCigar))
    
    def testCactusRealignSplitSequences(self):
        """Runs cactus realign, splitting indels longer than 100bp, and check
        that the coverage from the results is the same as the coverage from
        realigning with no arguments.."""
        for seqFile1, seqFile2 in seqFilePairGenerator():
            # Drop the lastz command since it's not needed. But this
            # is still convenient to use the same parameters as all
            # the other tests
            realignCommand, _ = getCommands(seqFile1, seqFile2)
            splitRealignCommand = realignCommand + " --splitIndelsLongerThanThis 100"
            realignOutput = getTempFile()
            splitRealignOutput = getTempFile()
            realignCommand += " > %s" % realignOutput
            splitRealignCommand += " > %s" % splitRealignOutput
            system(realignCommand)
            system(splitRealignCommand)
            # Check coverage on seqFile1
            splitRealignCoverage = popenCatch("cactus_coverage %s %s" % (seqFile1, splitRealignOutput))
            realignCoverage = popenCatch("cactus_coverage %s %s" % (seqFile1, realignOutput))
            self.assertTrue(splitRealignCoverage == realignCoverage)
            # Check coverage on seqFile2
            splitRealignCoverage = popenCatch("cactus_coverage %s %s" % (seqFile2, splitRealignOutput))
            realignCoverage = popenCatch("cactus_coverage %s %s" % (seqFile2, realignOutput))
            self.assertTrue(splitRealignCoverage == realignCoverage)
            os.remove(realignOutput)
            os.remove(splitRealignOutput)

    def testCactusRealignRescoreByIdentityAndProb(self):
        """Runs cactus realign using the default parameters and checks that the realigned output cigars align 
        the same subsequences.
        """
        for seqFile1, seqFile2 in seqFilePairGenerator():
            realignCommandByIdentity, lastzCommand = getCommands(seqFile1, seqFile2, realignArguments="--rescoreByIdentity")
            realignCommandByPosteriorProb = getCommands(seqFile1, seqFile2, realignArguments="--rescoreByPosteriorProb")[0]
            realignCommandByIdentityIgnoringGaps = getCommands(seqFile1, seqFile2, realignArguments="--rescoreByIdentityIgnoringGaps")[0]
            for realignLineByIdentity, realignLineByPosteriorProb, realignLineByIdentityIgnoringGaps, lastzLine in \
                                          zip([ i for i in popenCatch(realignCommandByIdentity).split("\n") if i != '' ], \
                                              [ i for i in popenCatch(realignCommandByPosteriorProb).split("\n") if i != '' ], \
                                              [ i for i in popenCatch(realignCommandByIdentityIgnoringGaps).split("\n") if i != '' ], \
                                              [ i for i in popenCatch(lastzCommand).split("\n") if i != '' ]):
                realignCigarByIdentity = cigarReadFromString(realignLineByIdentity)
                realignCigarByPosteriorProb = cigarReadFromString(realignLineByPosteriorProb)
                realignCigarByIdentityIgnoringGaps = cigarReadFromString(realignLineByIdentityIgnoringGaps)
                lastzCigar = cigarReadFromString(lastzLine)
                #Check scores are as expected
                self.assertTrue(realignCigarByIdentity.score >= 0)
                self.assertTrue(realignCigarByIdentity.score <= 100.0)
                self.assertTrue(realignCigarByPosteriorProb.score >= 0)
                self.assertTrue(realignCigarByPosteriorProb.score <= 100.0)
                self.assertTrue(realignCigarByIdentityIgnoringGaps.score >= 0)
                self.assertTrue(realignCigarByIdentityIgnoringGaps.score <= 100.0)
                #print "Scores", "Rescore by identity", realignCigarByIdentity.score, "Rescore by posterior prob", realignCigarByPosteriorProb.score, "Rescore by identity ignoring gaps", realignCigarByIdentityIgnoringGaps.score, "Lastz score", lastzCigar.score

def getCommands(seqFile1, seqFile2, realignArguments="", lastzArguments="--ambiguous=iupac"):  
    lastzCommand = "cactus_lastz --format=cigar %s %s[multiple][nameparse=darkspace] %s[nameparse=darkspace]" % (lastzArguments, seqFile1, seqFile2)
    realignCommand = "%s | cactus_realign %s %s %s" % (lastzCommand, realignArguments, seqFile1, seqFile2)
    return realignCommand, lastzCommand
                            
def seqFilePairGenerator():
     ##Get sequences
    encodePath = os.path.join(TestStatus.getPathToDataSets(), "MAY-2005")
    encodeRegions = [ "ENm00" + str(i) for i in xrange(1,2) ] #, 2) ] #Could go to six
    species = ("human", "mouse") #, "dog")#, "chimp") 
    #Other species to try "rat", "monodelphis", "macaque", "chimp"
    for encodeRegion in encodeRegions:
        regionPath = os.path.join(encodePath, encodeRegion)
        for i in xrange(len(species)):
            species1 = species[i]
            for species2 in species[i+1:]:
                seqFile1 = os.path.join(regionPath, "%s.%s.fa" % (species1, encodeRegion))
                seqFile2 = os.path.join(regionPath, "%s.%s.fa" % (species2, encodeRegion))
                yield seqFile1, seqFile2
        
if __name__ == '__main__':
    unittest.main()
