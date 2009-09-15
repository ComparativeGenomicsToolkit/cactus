import unittest

import cactus_setupTest
import cactus_alignerTest
import cactus_coreTest
import cactus_adjacencyTest
import cactus_workflowTest
import cactus_workflowCrashTest

from sonLib.bioio import parseSuiteTestOptions

def allSuites(): 
    cactus_setupSuite = unittest.makeSuite(cactus_setupTest.TestCase, 'test')
    cactus_alignerSuite = unittest.makeSuite(cactus_alignerTest.TestCase, 'test')
    cactus_coreSuite = unittest.makeSuite(cactus_coreTest.TestCase, 'test')
    cactus_adjacencySuite = unittest.makeSuite(cactus_adjacencyTest.TestCase, 'test')
    cactus_workflowSuite = unittest.makeSuite(cactus_workflowTest.TestCase, 'test')
    cactus_workflowCrashSuite = unittest.makeSuite(cactus_workflowCrashTest.TestCase, 'test')
    allTests = unittest.TestSuite((cactus_setupSuite, cactus_alignerSuite, 
                                   cactus_coreSuite, cactus_adjacencySuite, 
                                   cactus_workflowSuite, cactus_workflowCrashSuite))
    return allTests
        
def main():
    parseSuiteTestOptions()
    
    suite = allSuites()
    runner = unittest.TextTestRunner()
    runner.run(suite)
        
if __name__ == '__main__':
    main()