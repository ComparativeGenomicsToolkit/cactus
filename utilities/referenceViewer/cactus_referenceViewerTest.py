import unittest
import sys

from cactus.shared.test import parseCactusSuiteTestOptions
from sonLib.bioio import TestStatus

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import runWorkflow_multipleExamples

class TestCase(unittest.TestCase):
    def testCactus_Random(self):
        runWorkflow_multipleExamples(getCactusInputs_random, 
                                     testNumber=TestStatus.getTestSetup(), 
                                     buildReferencePDF=True)
        
def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()