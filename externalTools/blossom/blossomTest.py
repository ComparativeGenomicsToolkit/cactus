"""Tests the blossom program using blossomTest.py

Type python blossomTest.py --help to see options to the script.
"""

import unittest

import os
import sys

from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempDirectory
from sonLib.bioio import getTempFile
from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import saveInputs

class TestCase(unittest.TestCase):

    def setUp(self):
        #This is the number of random problems to solve, handed to the test code
        self.testNo = TestStatus.getTestSetup(shortTestNo=1, mediumTestNo=5, 
                                              longTestNo=10, veryLongTestNo=100)
        self.tempFiles = []
        self.tempDir = getTempDirectory(os.getcwd())
        self.tempBlossomDirectory = self.tempDir + "/tempBlossom"
        unittest.TestCase.setUp(self)

    def tearDown(self):
        for tempFile in self.tempFiles:
            os.remove(tempFile)
        unittest.TestCase.tearDown(self)
        # Comment out this line to avoid cleaning up all the files after the
        # test completes.
        system("rm -rf %s" % self.tempDir)

    def testBlossom(self):
        """ Tests blossom5 program using randGraph.py input
        """

        for test in xrange(self.testNo):
            tempInputFile = getTempFile()
            tempOutputFile = getTempFile()

            self.tempFiles.append(tempInputFile)
            self.tempFiles.append(tempOutputFile)

            # Create sample/test input graph file
            system("blossom_randGraph.py > %s" % tempInputFile)

            # Run blossom5
            system("blossom5 -e %s -w %s >& /dev/null" % (tempInputFile, tempOutputFile))

            # Now check if output is valid
            f = open(tempOutputFile, 'r')
            lineIdx = 0
            for line in f:
                line = line.rstrip()
                if lineIdx == 0:
                    (vertexNum, edgeNum) = line.split()
                    vertexNum = int(vertexNum)
                    edgeNum = int(edgeNum)
                    vertexArray = [0] * vertexNum

                    # Number of vertices must be even
                    self.assertEqual(vertexNum % 2, 0)

                    # Number of edges is half the number of vertices
                    self.assertEqual(vertexNum/2, edgeNum)
                else:
                    (vertexI, vertexJ,) = line.split()
                    vertexI = int(vertexI)
                    vertexJ = int(vertexJ)

                    vertexArray[vertexI] += 1
                    vertexArray[vertexJ] += 1

                    # Vertex indices must be 0<= i,j < V
                    self.assertTrue(vertexI in xrange(vertexNum))
                    self.assertTrue(vertexJ in xrange(vertexNum))
                lineIdx += 1

            # Must have the correct number of edges
            self.assertEqual(edgeNum, lineIdx-1)

            badCount = 0
            for i in vertexArray:
                if i != 1:
                    badCount += 1
            # Each vertex must be only in one edge
            self.assertEqual(badCount, 0)

            logger.info("Ran the test(s) of the blossom program okay")

def main():
    parseSuiteTestOptions()   
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()

