#!/usr/bin/env python

#Copyright (C) 2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest

from cactus.setup.cactus_setupTest import TestCase as setupTest
from cactus.blastAlignment.cactus_alignerTest import TestCase as alignerTest
from cactus.blastAlignment.cactus_batchTest import TestCase as batchTest
from cactus.pipeline.cactus_workflowTest import TestCase as workflowTest
from cactus.pipeline.cactus_evolverTest import TestCase as evolverTest
from cactus.baseAlignment.cactus_baseAlignerTest import TestCase as baseAlignerTest
from cactus.phylogeny.cactus_phylogenyTest import TestCase as phylogenyTest
from cactus.faces.cactus_fillAdjacenciesTest import TestCase as adjacenciesTest
from cactus.reference.cactus_referenceTest import TestCase as referenceTest
from cactus.hal.cactus_halTest import TestCase as halTest
from cactus.api.allTests import TestCase as aPITest
from cactus.externalTools.threeEdgeConnected.threeEdgeTests import TestCase as threeEdgeTest
from cactus.stPinchGraphs.allTests import TestCase as stPinchGraphsTest
from cactus.stCactusGraphs.allTests import TestCase as stCactusGraphsTest
from cactus.stCaf.allTests import TestCase as stCafTest
from cactus.normalisation.cactus_normalisationTest import TestCase as normalisationTest
from cactus.progressive.allTests import allSuites as progressiveSuite

from cactus.shared.test import parseCactusSuiteTestOptions

def allSuites(): 
    allTests = unittest.TestSuite((unittest.makeSuite(threeEdgeTest, 'test'),
                                   unittest.makeSuite(stPinchGraphsTest, 'test'),
                                   unittest.makeSuite(stCactusGraphsTest, 'test'),
                                   unittest.makeSuite(stCafTest, 'test'),
                                   unittest.makeSuite(setupTest, 'test'),
                                   unittest.makeSuite(alignerTest, 'test'),
                                   unittest.makeSuite(batchTest, 'test'),
                                   unittest.makeSuite(workflowTest, 'test'),
                                   unittest.makeSuite(evolverTest, 'test'),
                                   unittest.makeSuite(baseAlignerTest, 'test'),
                                   unittest.makeSuite(phylogenyTest, 'test'),
                                   unittest.makeSuite(adjacenciesTest, 'test'),
                                   unittest.makeSuite(referenceTest, 'test'),
                                   unittest.makeSuite(aPITest, 'test'),
                                   unittest.makeSuite(normalisationTest, 'test'),  
                                   unittest.makeSuite(halTest, 'test'),                                  
                                   progressiveSuite()))
    return allTests
        
def main():
    parseCactusSuiteTestOptions()
    
    suite = allSuites()
    runner = unittest.TextTestRunner()
    i = runner.run(suite)
    return len(i.failures) + len(i.errors)
        
if __name__ == '__main__':
    import sys
    sys.exit(main())
                
