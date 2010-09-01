import unittest
import sys

from cactus.shared.test import parseCactusSuiteTestOptions
from sonLib.bioio import TestStatus

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import runWorkflow_multipleExamples

class TestCase(unittest.TestCase):
    def testCactus_Random(self):
        runWorkflow_multipleExamples(getCactusInputs_random, 
                                     testNumber=TestStatus.getTestSetup(), 
                                     buildTrees=True, buildFaces=True, buildReference=False)
        
    def testCactus_Blanchette(self):
        runWorkflow_multipleExamples(getCactusInputs_blanchette,
                                     testRestrictions=(TestStatus.TEST_SHORT,), inverseTestRestrictions=True, 
                                     buildTrees=True, buildFaces=True, buildReference=False)
    
def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
