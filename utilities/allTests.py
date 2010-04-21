import unittest

from cactus.utilities.graphVizPlots.cactus_graphVizTest import TestCase as graphVizTest
from cactus.utilities.mAFs.cactus_MAFGeneratorTest import TestCase as mAFTest
from cactus.utilities.stats.cactus_treeStatsTest import TestCase as statsTest

from sonLib.bioio import parseSuiteTestOptions

def allSuites(): 
    allTests = unittest.TestSuite((unittest.makeSuite(graphVizTest, 'test'),
                                   unittest.makeSuite(mAFTest, 'test'),
                                   unittest.makeSuite(statsTest, 'test')))
    return allTests
        
def main():
    parseSuiteTestOptions()
     
    suite = allSuites()
    runner = unittest.TextTestRunner()
    runner.run(suite)
        
if __name__ == '__main__':
    main()
 