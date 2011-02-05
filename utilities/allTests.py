#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest

from cactus.utilities.graphVizPlots.cactus_graphVizTest import TestCase as graphVizTest
from cactus.utilities.mafs.cactus_MAFGeneratorTest import TestCase as mAFTest
from cactus.utilities.stats.cactus_treeStatsTest import TestCase as statsTest
from cactus.utilities.referenceViewer.cactus_referenceViewerTest import TestCase as referenceViewerTest

from sonLib.bioio import parseSuiteTestOptions

def allSuites(): 
    allTests = unittest.TestSuite((unittest.makeSuite(graphVizTest, 'test'),
                                   unittest.makeSuite(mAFTest, 'test'),
                                   unittest.makeSuite(statsTest, 'test'),
                                   unittest.makeSuite(referenceViewerTest, 'test')))
    return allTests
        
def main():
    parseSuiteTestOptions()
     
    suite = allSuites()
    runner = unittest.TextTestRunner()
    runner.run(suite)
        
if __name__ == '__main__':
    main()
 
