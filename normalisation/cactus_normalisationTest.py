import unittest
import sys

from sonLib.bioio import parseSuiteTestOptions
from sonLib.bioio import TestStatus

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import runWorkflow_multipleExamples

class TestCase(unittest.TestCase):
    def testCactusNormalisation_Random(self):
        runWorkflow_multipleExamples(getCactusInputs_random,
                                     testNumber=1000, #TestStatus.getTestSetup(),
                                     buildTrees=False, buildFaces=False, buildReference=False)
        
    def testCactusNormalisation_Blanchette(self):
        runWorkflow_multipleExamples(getCactusInputs_blanchette,
                                     testRestrictions=(TestStatus.TEST_SHORT,), inverseTestRestrictions=True,
                                     buildTrees=False, buildFaces=False, buildReference=False)
    
def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
