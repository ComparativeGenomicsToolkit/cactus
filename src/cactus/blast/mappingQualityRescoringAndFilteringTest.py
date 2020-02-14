import os
import pytest
import random
import time
import unittest

from toil.job import Job
from toil.common import Toil

from sonLib.bioio import getTempFile, getTempDirectory, system
from sonLib.bioio import TestStatus
from textwrap import dedent
from cactus.shared.common import cactus_call, runSelfLastz
from cactus.shared.test import getCactusInputs_encode
from cactus.blast.mappingQualityRescoringAndFiltering import mappingQualityRescoring
from cactus.shared.common import makeURL

from cactus.shared.test import getCactusInputs_evolverMammals
from cactus.shared.test import getCactusInputs_evolverPrimates
from setuptools.dist import sequence

@pytest.mark.blast
@TestStatus.needsTestData
class TestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        # simple test data
        self.simpleInputCigarPath = getTempFile()
        self.inputCigars = [
            "cigar: simpleSeqBC 0 9 + simpleSeqAC 10 0 - 5 M 8 D 1 M 1" ,
            "cigar: simpleSeqB1 18 9 - simpleSeqA1 6 2 - 4 M 3 I 5 M 1",
            "cigar: simpleSeqB1 32 30 - simpleSeqA2 7 9 + 72 M 2",
            "cigar: simpleSeqB1 9 18 + simpleSeqA1 2 6 + 1 M 3 I 5 M 1",
            "cigar: simpleSeqB1 18 28 + simpleSeqA2 0 10 + 8 M 1 I 2 M 2 D 2 M 5",
            "cigar: simpleSeqB1 28 30 + simpleSeqA2 6 8 + 3 M 2",
            "cigar: simpleSeqZ1 0 1 + simpleSeqA1 6 7 + 3 M 1",
            "cigar: simpleSeqC1 0 5 + simpleSeqD 0 5 + 2 M 5",
            "cigar: simpleSeqD 5 10 + simpleSeqC1 5 10 + 8 M 5",
            "cigar: simpleSeqC1 10 15 + simpleSeqC1 15 20 + 19 M 5",
            "cigar: simpleSeqNonExistent 0 10 + simpleSeqC1 0 10 + 0.5 M 10" ]

        with open(self.simpleInputCigarPath, 'w') as fH:
            fH.write("\n".join(self.inputCigars) + "\n")

        self.sortedNonOverlappingInputCigars = [
            'cigar: simpleSeqB1 9 18 + simpleSeqA1 2 6 + 1.000000 M 3 I 5 M 1',
            'cigar: simpleSeqB1 9 18 + simpleSeqA1 2 6 + 4.000000 M 1 I 5 M 3',
            'cigar: simpleSeqZ1 0 1 + simpleSeqA1 6 7 + 3.000000 M 1',
            'cigar: simpleSeqB1 18 24 + simpleSeqA2 0 6 + 8.000000 M 1 I 2 M 2 D 2 M 1',
            'cigar: simpleSeqB1 28 29 + simpleSeqA2 6 7 + 3.000000 M 1',
            'cigar: simpleSeqB1 24 25 + simpleSeqA2 6 7 + 8.000000 M 1',
            'cigar: simpleSeqB1 29 30 + simpleSeqA2 7 8 + 3.000000 M 1',
            'cigar: simpleSeqB1 32 31 - simpleSeqA2 7 8 + 72.000000 M 1',
            'cigar: simpleSeqB1 25 26 + simpleSeqA2 7 8 + 8.000000 M 1',
            'cigar: simpleSeqB1 31 30 - simpleSeqA2 8 9 + 72.000000 M 1',
            'cigar: simpleSeqB1 26 27 + simpleSeqA2 8 9 + 8.000000 M 1',
            'cigar: simpleSeqB1 27 28 + simpleSeqA2 9 10 + 8.000000 M 1',
            'cigar: simpleSeqBC 9 0 - simpleSeqAC 0 10 + 5.000000 M 1 D 1 M 8',
            'cigar: simpleSeqA1 2 6 + simpleSeqB1 9 18 + 4.000000 M 1 D 5 M 3',
            'cigar: simpleSeqA1 2 6 + simpleSeqB1 9 18 + 1.000000 M 3 D 5 M 1',
            'cigar: simpleSeqA2 0 10 + simpleSeqB1 18 28 + 8.000000 M 1 D 2 M 2 I 2 M 5',
            'cigar: simpleSeqA2 6 8 + simpleSeqB1 28 30 + 3.000000 M 2',
            'cigar: simpleSeqA2 9 7 - simpleSeqB1 30 32 + 72.000000 M 2',
            'cigar: simpleSeqAC 10 0 - simpleSeqBC 0 9 + 5.000000 M 8 I 1 M 1',
            'cigar: simpleSeqD 0 5 + simpleSeqC1 0 5 + 2.000000 M 5',
            'cigar: simpleSeqNonExistent 0 5 + simpleSeqC1 0 5 + 0.500000 M 5',
            'cigar: simpleSeqNonExistent 5 10 + simpleSeqC1 5 10 + 0.500000 M 5',
            'cigar: simpleSeqD 5 10 + simpleSeqC1 5 10 + 8.000000 M 5',
            'cigar: simpleSeqC1 15 20 + simpleSeqC1 10 15 + 19.000000 M 5',
            'cigar: simpleSeqC1 10 15 + simpleSeqC1 15 20 + 19.000000 M 5',
            'cigar: simpleSeqC1 0 5 + simpleSeqD 0 5 + 2.000000 M 5',
            'cigar: simpleSeqC1 5 10 + simpleSeqD 5 10 + 8.000000 M 5',
            'cigar: simpleSeqC1 0 10 + simpleSeqNonExistent 0 10 + 0.500000 M 10',
            'cigar: simpleSeqA1 6 7 + simpleSeqZ1 0 1 + 3.000000 M 1',
        ]

        self.filteredSortedNonOverlappingInputCigars = [
            'cigar: simpleSeqB1 9 18 + simpleSeqA1 2 6 + 30.004341 M 1 I 5 M 3',
            'cigar: simpleSeqZ1 0 1 + simpleSeqA1 6 7 + 60.000000 M 1',
            'cigar: simpleSeqB1 18 24 + simpleSeqA2 0 6 + 60.000000 M 1 I 2 M 2 D 2 M 1',
            'cigar: simpleSeqB1 24 25 + simpleSeqA2 6 7 + 50.000042 M 1',
            'cigar: simpleSeqB1 32 31 - simpleSeqA2 7 8 + 60.000000 M 1',
            'cigar: simpleSeqB1 31 30 - simpleSeqA2 8 9 + 60.000000 M 1',
            'cigar: simpleSeqB1 27 28 + simpleSeqA2 9 10 + 60.000000 M 1',
            'cigar: simpleSeqBC 9 0 - simpleSeqAC 0 10 + 60.000000 M 1 D 1 M 8',
            'cigar: simpleSeqA1 2 6 + simpleSeqB1 9 18 + 30.004341 M 1 D 5 M 3',
            'cigar: simpleSeqA2 0 10 + simpleSeqB1 18 28 + 60.000000 M 1 D 2 M 2 I 2 M 5',
            'cigar: simpleSeqA2 6 8 + simpleSeqB1 28 30 + 60.000000 M 2',
            'cigar: simpleSeqA2 9 7 - simpleSeqB1 30 32 + 60.000000 M 2',
            'cigar: simpleSeqAC 10 0 - simpleSeqBC 0 9 + 60.000000 M 8 I 1 M 1',
            'cigar: simpleSeqD 0 5 + simpleSeqC1 0 5 + 15.135209 M 5',
            'cigar: simpleSeqD 5 10 + simpleSeqC1 5 10 + 60.000000 M 5',
            'cigar: simpleSeqC1 15 20 + simpleSeqC1 10 15 + 60.000000 M 5',
            'cigar: simpleSeqC1 10 15 + simpleSeqC1 15 20 + 60.000000 M 5',
            'cigar: simpleSeqC1 0 5 + simpleSeqD 0 5 + 60.000000 M 5',
            'cigar: simpleSeqC1 5 10 + simpleSeqD 5 10 + 60.000000 M 5',
            'cigar: simpleSeqC1 0 10 + simpleSeqNonExistent 0 10 + 60.000000 M 10',
            'cigar: simpleSeqA1 6 7 + simpleSeqZ1 0 1 + 60.000000 M 1'
        ]

        self.simpleOutputCigarPath = getTempFile()
        self.simpleOutputCigarPath2 = getTempFile()
        self.sequenceFilePath = getTempFile()
        self.logLevelString = "DEBUG"

        self.tempDir = getTempDirectory(os.getcwd())

    def tearDown(self):
        unittest.TestCase.tearDown(self)
        os.remove(self.simpleInputCigarPath)
        os.remove(self.simpleOutputCigarPath)
        os.remove(self.simpleOutputCigarPath2)
        os.remove(self.sequenceFilePath)
        system("rm -rf %s" % self.tempDir)

    @staticmethod
    def makeCigar(coordinates1, coordinates2, score, ops):
        cigar = "cigar: %s %s %f %s" % (" ".join(map(str, coordinates2)),
                                        " ".join(map(str, coordinates1)),
                                        float(score), " ".join(map(str, ops)))
        return cigar

    @TestStatus.shortLength
    def testMirrorAndOrientAlignments(self):
        cactus_call(parameters=["cactus_mirrorAndOrientAlignments",
                                 self.logLevelString,
                                 self.simpleInputCigarPath,
                                 self.simpleOutputCigarPath])
        with open(self.simpleOutputCigarPath, 'r') as fh:
            outputCigars = [ cigar[:-1] for cigar in fh.readlines() ] # Remove new lines

        # For each input alignment check that we have the two, oriented alignments
        for inputCigar in self.inputCigars:
            name1, start1, end1, strand1 = inputCigar.split()[5:9]
            start1, end1 = int(start1), int(end1)
            coordinates1 = name1, start1, end1, strand1

            name2, start2, end2, strand2 = inputCigar.split()[1:5]
            start2, end2 = int(start2), int(end2)
            coordinates2 = name2, start2, end2, strand2

            score = inputCigar.split()[9]
            ops = inputCigar.split()[10:]

            def invertStrand(coordinates):
                # cigar: simpleSeqB1 0 9 + simpleSeqA1 10 0 - 0 M 8 D 1 M 1
                # cigar: simpleSeqB1 9 0 + simpleSeqA1 0 10 - 0 M 1 D 1 M 8
                name, start, end, strand = coordinates
                assert strand in ("+", "-")
                if strand == "+":
                    return name, end, start, "-"
                return name, end, start, "+"

            def reverseOps(ops):
                l = ops[:]
                l.reverse()
                l2 = []
                for i, j in zip(l[1::2], l[::2]):
                    l2 += [ i, j ]
                return l2

            def invertOpStrands(ops):
                l = [ "I" if op == "D" else ("D" if op == "I" else op) for op in ops[::2] ]
                l2 = []
                for op, length in zip(l, ops[1::2]):
                    l2 += [ op, length ]
                return l2

            if strand1 == "+":
                self.assertTrue(self.makeCigar(coordinates1, coordinates2, score, ops) in outputCigars)
            else:
                # Invert the strands
                self.assertTrue(self.makeCigar(invertStrand(coordinates1),
                                          invertStrand(coordinates2), score, reverseOps(ops)) in outputCigars)

            if strand2 == "+":
                self.assertTrue(self.makeCigar(coordinates2, coordinates1, score,
                                          invertOpStrands(ops)) in outputCigars)
            else:
                self.assertTrue(self.makeCigar(invertStrand(coordinates2), invertStrand(coordinates1),
                                            score, invertOpStrands(reverseOps(ops))) in outputCigars)

    @TestStatus.shortLength
    def testSplitAlignmentsOverlaps(self):
        self.inputCigars = [
            'cigar: simpleSeqB1 9 18 + simpleSeqA1 2 6 + 1.000000 M 3 I 5 M 1',
            'cigar: simpleSeqB1 9 18 + simpleSeqA1 2 6 + 4.000000 M 1 I 5 M 3',
            'cigar: simpleSeqZ1 0 1 + simpleSeqA1 6 7 + 3.000000 M 1',
            'cigar: simpleSeqB1 18 28 + simpleSeqA2 0 10 + 8.000000 M 1 I 2 M 2 D 2 M 5',
            'cigar: simpleSeqB1 28 30 + simpleSeqA2 6 8 + 3.000000 M 2',
            'cigar: simpleSeqB1 32 30 - simpleSeqA2 7 9 + 72.000000 M 2',
            'cigar: simpleSeqBC 9 0 - simpleSeqAC 0 10 + 5.000000 M 1 D 1 M 8',
            'cigar: simpleSeqA1 2 6 + simpleSeqB1 9 18 + 1.000000 M 3 D 5 M 1',
            'cigar: simpleSeqA1 2 6 + simpleSeqB1 9 18 + 4.000000 M 1 D 5 M 3',
            'cigar: simpleSeqA2 0 10 + simpleSeqB1 18 28 + 8.000000 M 1 D 2 M 2 I 2 M 5',
            'cigar: simpleSeqA2 6 8 + simpleSeqB1 28 30 + 3.000000 M 2',
            'cigar: simpleSeqA2 9 7 - simpleSeqB1 30 32 + 72.000000 M 2',
            'cigar: simpleSeqAC 10 0 - simpleSeqBC 0 9 + 5.000000 M 8 I 1 M 1',
            'cigar: simpleSeqD 0 5 + simpleSeqC1 0 5 + 2.000000 M 5',
            'cigar: simpleSeqNonExistent 0 10 + simpleSeqC1 0 10 + 0.500000 M 10',
            'cigar: simpleSeqD 5 10 + simpleSeqC1 5 10 + 8.000000 M 5',
            'cigar: simpleSeqC1 15 20 + simpleSeqC1 10 15 + 19.000000 M 5',
            'cigar: simpleSeqC1 10 15 + simpleSeqC1 15 20 + 19.000000 M 5',
            'cigar: simpleSeqC1 0 5 + simpleSeqD 0 5 + 2.000000 M 5',
            'cigar: simpleSeqC1 5 10 + simpleSeqD 5 10 + 8.000000 M 5',
            'cigar: simpleSeqC1 0 10 + simpleSeqNonExistent 0 10 + 0.500000 M 10',
            'cigar: simpleSeqA1 6 7 + simpleSeqZ1 0 1 + 3.000000 M 1'
        ]
        with open(self.simpleInputCigarPath, 'w') as fH:
            fH.write("\n".join(self.inputCigars) + "\n")

        cactus_call(parameters=["cactus_splitAlignmentOverlaps",
                                 self.logLevelString,
                                 self.simpleInputCigarPath,
                                 self.simpleOutputCigarPath])

        with open(self.simpleOutputCigarPath, 'r') as fh:
            outputCigars = [ cigar[:-1] for cigar in fh.readlines() ] # Remove new lines

        # Get start and end coordinates of cigars
        ends = set()
        for inputCigar in self.inputCigars:
            name1, start1, end1, strand1 = inputCigar.split()[5:9]
            ends.add((name1, int(start1)))
            ends.add((name1, int(end1)))
            assert strand1 == "+"

        # Count of expected number of chopped up cigars
        totalExpectedCigars = 0

        # Function to split a list of ops into a prefix and suffix list
        def splitPrefixOps(ops, cutPoint):
            pOps, sOps = [], []
            j = 0
            for i in range(0, len(ops), 2):
                op, length = ops[i], int(ops[i+1])
                assert op in ("I", "D", "M")
                if op == "I":
                    pOps.append(op)
                    pOps.append(length)
                    continue
                if j + length <= cutPoint:
                    pOps.append(op)
                    pOps.append(length)
                    j += length
                    if j == cutPoint:
                        break
                else:
                    assert j + length > cutPoint
                    pOps.append(op)
                    pOps.append(cutPoint - j)
                    sOps.append(op)
                    sOps.append(length - (cutPoint - j))
                    break
            sOps += ops[i+2:]

            return pOps, sOps

        # For each cigar:
        for inputCigar in self.inputCigars:
            name1, start1, end1, strand1 = inputCigar.split()[5:9]
            start1, end1 = int(start1), int(end1)
            assert strand1 == "+"
            name2, start2, end2, strand2 = inputCigar.split()[1:5]
            start2, end2 = int(start2), int(end2)
            score = float(inputCigar.split()[9])
            ops = inputCigar.split()[10:]

            # For each intermediate chop point
            i = start1
            for j in range(start1+1, end1+1):
                if (name1, j) in ends:
                    # Chop up cigar
                    coordinates1 = name1, i, j, "+"

                    # Get sublist of ops
                    pOps, subOps = splitPrefixOps(ops, i - start1)
                    subOps, sOps = splitPrefixOps(subOps, j - i)

                    x = lambda ops : sum([ int(ops[k+1]) for k in range(0, len(ops), 2) if ops[k] != 'D' ])
                    k = x(pOps)
                    l = k + x(subOps)

                    # Get second coordinates
                    if strand2 == "+":
                        coordinates2 = name2, start2 + k, start2 + l, strand2
                    else:
                        assert strand2 == "-"
                        coordinates2 = name2, start2 - k, start2 - l, strand2

                    choppedCigar = self.makeCigar(coordinates1, coordinates2, score, subOps)

                    # Check each chopped up cigar is in output
                    self.assertTrue(choppedCigar in outputCigars)

                    # Inc. number of expected cigars
                    totalExpectedCigars += 1

                    # Check previous coordinate
                    i = j

        # Check we have the expected number of cigars
        self.assertEqual(totalExpectedCigars, len(outputCigars))

    @TestStatus.shortLength
    def testCalculateMappingQualities(self):
        with open(self.simpleInputCigarPath, 'w') as fH:
            fH.write("\n".join(self.sortedNonOverlappingInputCigars) + "\n")

        cactus_call(parameters=[ "cactus_calculateMappingQualities",
                                 self.logLevelString,
                                 '1', '0', "1.0",
                                 self.simpleOutputCigarPath,
                                 self.simpleInputCigarPath ])

        with open(self.simpleOutputCigarPath, 'r') as fh:
            outputCigars = [ cigar[:-1] for cigar in fh.readlines() ] # Remove new lines

        self.assertEqual(self.filteredSortedNonOverlappingInputCigars, outputCigars)

    def runToilPipeline(self, alignmentsFile, alpha=0.001):
        # Tests the toil pipeline
        options = Job.Runner.getDefaultOptions(os.path.join(self.tempDir, "toil"))
        options.logLevel = self.logLevelString

        with Toil(options) as toil:
            # Import the input file into the job store
            inputAlignmentFileID = toil.importFile(makeURL(alignmentsFile))

            rootJob = Job.wrapJobFn(mappingQualityRescoring, inputAlignmentFileID,
                                    minimumMapQValue=0, maxAlignmentsPerSite=1, alpha=alpha, logLevel=self.logLevelString)

            primaryOutputAlignmentsFileID, secondaryOutputAlignmentsFileID = toil.start(rootJob)
            toil.exportFile(primaryOutputAlignmentsFileID, makeURL(self.simpleOutputCigarPath))
            toil.exportFile(secondaryOutputAlignmentsFileID, makeURL(self.simpleOutputCigarPath2))

        # Check output
        with open(self.simpleOutputCigarPath, 'r') as fh:
            primaryOutputCigars = [ cigar[:-1] for cigar in fh.readlines() ] # Remove new lines

        with open(self.simpleOutputCigarPath2, 'r') as fh:
            secondaryOutputCigars = [ cigar[:-1] for cigar in fh.readlines() ] # Remove new lines

        return primaryOutputCigars + secondaryOutputCigars

    @TestStatus.shortLength
    def testMappingQualityRescoringAndFiltering(self):
        """
        Tests the complete pipeline.
        """
        outputCigars = self.runToilPipeline(self.simpleInputCigarPath, alpha=1.0)

        self.assertEqual(self.filteredSortedNonOverlappingInputCigars, outputCigars)

    def alignAndRunPipeline(self, concatenatedSequenceFile):
        # Run lastz
        startTime = time.time()
        runSelfLastz(concatenatedSequenceFile, self.simpleInputCigarPath, "")
        print(("It took %s seconds to run lastz" % (time.time() - startTime)))
        with open(self.simpleInputCigarPath, 'r') as fh:
            inputCigars = [ cigar[:-1] for cigar in fh.readlines() ] # Remove new lines
        print(("There are %s cigars from lastz" % len(inputCigars)))

        # Run toil pipeline
        startTime = time.time()
        outputCigars = self.runToilPipeline(self.simpleInputCigarPath)
        print(("It took %s seconds to run unique mapping pipeline" % (time.time() - startTime)))

        # Check output
        print(("Total input cigars:", len(inputCigars), "total output cigars", len(outputCigars)))
        #print outputCigars

    @TestStatus.longLength
    def testMappingQualityRescoringAndFilteringScaled_Mammals(self):
        # Test with a larger number of realistic input cigars

        # Get inputs
        sequenceFiles, newickTreeString = getCactusInputs_evolverMammals()
        concatenatedSequenceFile = os.path.join(self.tempDir, "tempMammals.fa")
         # Cat the sequences together
        system('cat %s > %s' % (" ".join(sequenceFiles), concatenatedSequenceFile))

        # Align and run the pipeline
        self.alignAndRunPipeline(concatenatedSequenceFile)

        os.remove(concatenatedSequenceFile) # This is a hack, because the local mode of lastz call is not working

    @TestStatus.mediumLength
    def testMappingQualityRescoringAndFilteringScaled_Primates(self):
        # Test with a larger number of realistic input cigars

        # Get inputs
        sequenceFiles, newickTreeString = getCactusInputs_evolverPrimates()
        concatenatedSequenceFile = "./tempPrimates.fa" #os.path.join(self.tempDir, "tempPrimates.fa")
         # Cat the sequences together
        system('cat %s > %s' % (" ".join(sequenceFiles), concatenatedSequenceFile))

        # Align and run the pipeline
        self.alignAndRunPipeline(concatenatedSequenceFile)

        os.remove(concatenatedSequenceFile) # This is a hack, because the local mode of lastz call is not working

if __name__ == '__main__':
    unittest.main()
