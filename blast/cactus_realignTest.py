import unittest

from sonLib.bioio import cigarRead
from sonLib.bioio import system
from sonLib.bioio import logger

class TestCase(unittest.TestCase):
    
    def setUp(self):
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        
    def testCactusRealignDummy(self):
        """Runs cactus realign using the "dummy" mode
        and checks the output is equivalent to what you'd get by just running lastz.
        """
        for seqFile1, seqFile2 in seqFilePairGenerator():
            realignCommand, lastzCommand = getCommands(seqFile1, seqFile2, "--dummy")
            for realignLine, lastzLine in zip([ i for i in popenCatch(realignCommand).split("\n") if i != '' ], 
                                              [ i for i in popenCatch(lastzCommand).split("\n") if i != '' ]):
                realignCigar = cigarRead(realignLine)
                lastzCigar = cigarRead(lastzLine)
                self.assertTrue(realignCigar == lastzCigar)
    
    def testCactusRealign(self):
        """Runs cactus realign using the default parameters and checks that the realigned output cigars align 
        the same subsequences.
        """
        for seqFile1, seqFile2 in seqFilePairGenerator():
            realignCommand, lastzCommand = getCommands(seqFile1, seqFile2)
            for realignLine, lastzLine in zip([ i for i in popenCatch(realignCommand).split("\n") if i != '' ], 
                                              [ i for i in popenCatch(lastzCommand).split("\n") if i != '' ]):
                realignCigar = cigarRead(realignLine)
                lastzCigar = cigarRead(lastzLine)
                self.assertTrue(realignCigar.sameCoordinates(lastzCigar))

def getCommands(seqFile1, seqFile2, realignArguments="", lastzArguments="--hspthresh=1800 --ambiguous=iupac"):  
    lastzArguments="--hspthresh=1800 --ambiguous=iupac"
    realignArguments="--dummy"
    lastzCommand = "lastz --format=cigar %s %s[multiple][nameparse=darkspace] %s[nameparse=darkspace]" % (lastzArguments, seqFile1, seqFile2)
    realignCommand = "%s | cactus_realign %s %s %s" % (lastzCommand, realignArguments, seqFile1, seqFile2)
    return realignCommand, lastzCommand
                            
def seqPairGenerator():
     ##Get sequences
    encodePath = os.path.join(TestStatus.getPathToDataSets(), "MAY-2005")
    encodeRegions = [ "ENm00" + str(i) for i in xrange(1,2) ] #, 2) ] #Could go to six
    species = ("human", "mouse", "dog")
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