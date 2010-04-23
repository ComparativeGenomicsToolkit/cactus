import unittest
import sys
import os

from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import system
from sonLib.bioio import getTempDirectory

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import runWorkflow_multipleExamples

class TestCase(unittest.TestCase):
    def testCactus_Random(self):
        tempOutputDir = getTempDirectory(os.getcwd())
        runWorkflow_multipleExamples(getCactusInputs_random, 
                                     outputDir=tempOutputDir,
                                     testNumber=TestStatus.getTestSetup(), 
                                     makeCactusTreeStats=True)
        #Now check we can build a latex stats file..
        l = []
        for test in xrange(TestStatus.getTestSetup()):
            l.append(os.path.join(tempOutputDir, str(test), "cactusStats.xml"))
            l.append("region%s" % test)
        statsFileTEX = os.path.join(tempOutputDir, "cactusStats.tex")
        system("cactus_treeStatsToLatexTables.py --outputFile %s %s" % \
                (statsFileTEX, " ".join(l)))
        system("rm -rf %s" % tempOutputDir)
        
def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()