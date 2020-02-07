import os
import pytest
import sys
import unittest

from sonLib.bioio import cigarReadFromString, cigarWrite
from sonLib.bioio import popenCatch
from sonLib.bioio import logger
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempFile
from sonLib.bioio import getTempDirectory
from sonLib.bioio import system

from cactus.shared.common import runLastz
from cactus.shared.common import runSelfLastz
from cactus.shared.common import runCactusRealign
from cactus.shared.common import runCactusSelfRealign
from cactus.shared.common import runCactusCoverage

@pytest.mark.blast
@TestStatus.needsTestData
class TestCase(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)
        self.tempDir = getTempDirectory(os.getcwd())
        self.tempFiles = []
        unittest.TestCase.setUp(self)
        self.tempOutputFile = os.path.join(self.tempDir, "results1.txt")
        self.tempFiles.append(self.tempOutputFile)
        self.tempOutputFile2 = os.path.join(self.tempDir, "results2.txt")
        self.tempFiles.append(self.tempOutputFile2)
        self.encodePath = os.path.join(TestStatus.getPathToDataSets(), "MAY-2005")
        self.defaultLastzArguments = "--ambiguous=iupac"
        self.defaultRealignArguments = ""

    def tearDown(self):
        unittest.TestCase.tearDown(self)

    @TestStatus.mediumLength
    def testCactusRealignDummy(self):
        """Runs cactus realign using the "rescoreOriginalAlignment" mode
        and checks the output is equivalent to what you'd get by just running lastz.
        """
        for seqFile1, seqFile2 in seqFilePairGenerator():

            lastzOutput = getTempFile(rootDir=self.tempDir)
            runLastz(seqFile1, seqFile2, alignmentsFile=lastzOutput,
                     lastzArguments=self.defaultLastzArguments)
            realignOutput = getTempFile(rootDir=self.tempDir)
            runCactusRealign(seqFile1, seqFile2, inputAlignmentsFile = lastzOutput,
                             outputAlignmentsFile = realignOutput,
                             realignArguments=self.defaultRealignArguments + " --rescoreOriginalAlignment")

            for realignLine, lastzLine in zip([ i for i in open(lastzOutput, 'r') if i != '' ],
                                              [ i for i in open(realignOutput, 'r') if i != '' ]):
                realignCigar = cigarReadFromString(realignLine)
                lastzCigar = cigarReadFromString(lastzLine)
                self.assertTrue(realignCigar != None)
                self.assertTrue(realignCigar == lastzCigar)

    @TestStatus.mediumLength
    def testCactusRealign(self):
        """Runs cactus realign using the default parameters and checks that the realigned output cigars align
        the same subsequences.
        """
        for seqFile1, seqFile2 in seqFilePairGenerator():
            lastzOutput = getTempFile(rootDir=self.tempDir)
            runLastz(seqFile1, seqFile2, alignmentsFile=lastzOutput,
                     lastzArguments=self.defaultLastzArguments)
            realignOutput = getTempFile(rootDir=self.tempDir)
            runCactusRealign(seqFile1, seqFile2, inputAlignmentsFile = lastzOutput,
                             outputAlignmentsFile = realignOutput,
                             realignArguments=self.defaultRealignArguments)

            for realignLine, lastzLine in zip([ i for i in open(lastzOutput, 'r') if i != '' ],
                                              [ i for i in open(realignOutput, 'r') if i != '' ]):
                realignCigar = cigarReadFromString(realignLine)
                lastzCigar = cigarReadFromString(lastzLine)
                self.assertTrue(realignCigar.sameCoordinates(lastzCigar))

    @TestStatus.mediumLength
    def testCactusRealignSplitSequences(self):
        """Runs cactus realign, splitting indels longer than 100bp, and check
        that the coverage from the results is the same as the coverage from
        realigning with no arguments.."""
        for seqFile1, seqFile2 in seqFilePairGenerator():
            lastzOutput = getTempFile(rootDir=self.tempDir)
            runLastz(seqFile1, seqFile2, alignmentsFile=lastzOutput,
                     lastzArguments=self.defaultLastzArguments)

            realignOutput = getTempFile(rootDir=self.tempDir)
            runCactusRealign(seqFile1, seqFile2, inputAlignmentsFile=lastzOutput,
                             outputAlignmentsFile=realignOutput,
                             realignArguments=self.defaultRealignArguments)

            splitRealignOutput = getTempFile(rootDir=self.tempDir)
            runCactusRealign(seqFile1, seqFile2, inputAlignmentsFile=lastzOutput,
                             outputAlignmentsFile=splitRealignOutput,
                             realignArguments=self.defaultRealignArguments + " --splitIndelsLongerThanThis 100")

            # Check coverage on seqFile1
            splitRealignCoverage = runCactusCoverage(seqFile1, splitRealignOutput)
            realignCoverage = runCactusCoverage(seqFile1, realignOutput)
            self.assertTrue(splitRealignCoverage == realignCoverage)
            # Check coverage on seqFile2
            splitRealignCoverage = runCactusCoverage(seqFile2, splitRealignOutput)
            realignCoverage = runCactusCoverage(seqFile2, realignOutput)
            self.assertTrue(splitRealignCoverage == realignCoverage)
            os.remove(realignOutput)
            os.remove(splitRealignOutput)

    @TestStatus.mediumLength
    def testCactusRealignRescoreByIdentityAndProb(self):
        """Runs cactus realign using the default parameters and checks that the realigned output cigars align
        the same subsequences.
        """
        for seqFile1, seqFile2 in seqFilePairGenerator():
            lastzOutput = getTempFile(rootDir=self.tempDir)
            runLastz(seqFile1, seqFile2, alignmentsFile=lastzOutput,
                     lastzArguments=self.defaultLastzArguments)

            realignByIdentityOutput = getTempFile(rootDir=self.tempDir)
            runCactusRealign(seqFile1, seqFile2, inputAlignmentsFile=lastzOutput,
                             outputAlignmentsFile=realignByIdentityOutput,
                             realignArguments=self.defaultRealignArguments + " --rescoreByIdentity")

            realignByPosteriorProbOutput = getTempFile(rootDir=self.tempDir)
            runCactusRealign(seqFile1, seqFile2, inputAlignmentsFile=lastzOutput,
                             outputAlignmentsFile=realignByPosteriorProbOutput,
                             realignArguments=self.defaultRealignArguments + " --rescoreByPosteriorProb")

            realignByIdentityIgnoringGapsOutput = getTempFile(rootDir=self.tempDir)
            runCactusRealign(seqFile1, seqFile2, inputAlignmentsFile=lastzOutput,
                             outputAlignmentsFile=realignByIdentityIgnoringGapsOutput,
                             realignArguments=self.defaultRealignArguments + " --rescoreByIdentityIgnoringGaps")
            for realignLineByIdentity, realignLineByPosteriorProb, realignLineByIdentityIgnoringGaps, lastzLine in \
                                          zip([ i for i in open(realignByIdentityOutput, 'r') if i != '' ], \
                                              [ i for i in open(realignByPosteriorProbOutput, 'r') if i != '' ], \
                                              [ i for i in open(realignByIdentityIgnoringGapsOutput, 'r') if i != '' ], \
                                              [ i for i in open(lastzOutput, 'r') if i != '' ]):
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


def seqFilePairGenerator():
     ##Get sequences
    encodePath = os.path.join(TestStatus.getPathToDataSets(), "MAY-2005")
    encodeRegions = [ "ENm00" + str(i) for i in range(1,2) ] #, 2) ] #Could go to six
    species = ("human", "mouse") #, "dog")#, "chimp")
    #Other species to try "rat", "monodelphis", "macaque", "chimp"
    for encodeRegion in encodeRegions:
        regionPath = os.path.join(encodePath, encodeRegion)
        for i in range(len(species)):
            species1 = species[i]
            for species2 in species[i+1:]:
                seqFile1 = os.path.join(regionPath, "%s.%s.fa" % (species1, encodeRegion))
                seqFile2 = os.path.join(regionPath, "%s.%s.fa" % (species2, encodeRegion))
                yield seqFile1, seqFile2

if __name__ == '__main__':
    unittest.main()
