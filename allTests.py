#!/usr/bin/env python

#Copyright (C) 2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest

from cactus.setup.cactus_setupTest import TestCase as setupTest
from cactus.blastAlignment.cactus_alignerTest import TestCase as alignerTest
from cactus.blastAlignment.cactus_batchTest import TestCase as batchTest
from cactus.core.cactus_coreTest import TestCase as coreTest
from cactus.pipeline.cactus_workflowTest import TestCase as workflowTest
from cactus.pipeline.cactus_evolverTest import TestCase as evolverTest
from cactus.baseAlignment.cactus_baseAlignerTest import TestCase as baseAlignerTest
from cactus.phylogeny.cactus_phylogenyTest import TestCase as phylogenyTest
from cactus.faces.cactus_fillAdjacenciesTest import TestCase as adjacenciesTest
from cactus.reference.cactus_referenceTest import TestCase as referenceTest
from cactus.api.allTests import TestCase as aPITest
from cactus.externalTools.threeEdgeConnected.threeEdgeTests import TestCase as threeEdgeTest
from cactus.externalTools.matchGraph.matchGraphTest import TestCase as matchGraphTest
from cactus.externalTools.blossom.blossomTest import TestCase as blossomTest
from cactus.normalisation.cactus_normalisationTest import TestCase as normalisationTest
 
from cactus.shared.test import parseCactusSuiteTestOptions

def allSuites(): 
    allTests = unittest.TestSuite((unittest.makeSuite(threeEdgeTest, 'test'),
                                   unittest.makeSuite(matchGraphTest, 'test'),
                                   unittest.makeSuite(blossomTest, 'test'),
                                   unittest.makeSuite(setupTest, 'test'),
                                   unittest.makeSuite(alignerTest, 'test'),
                                   unittest.makeSuite(batchTest, 'test'),
                                   unittest.makeSuite(coreTest, 'test'),
                                   unittest.makeSuite(workflowTest, 'test'),
                                   unittest.makeSuite(evolverTest, 'test'),
                                   unittest.makeSuite(baseAlignerTest, 'test'),
                                   unittest.makeSuite(phylogenyTest, 'test'),
                                   unittest.makeSuite(adjacenciesTest, 'test'),
                                   unittest.makeSuite(referenceTest, 'test'),
                                   unittest.makeSuite(aPITest, 'test'),
                                   unittest.makeSuite(normalisationTest, 'test')))
    return allTests
        
def main():
    parseCactusSuiteTestOptions()
    
    suite = allSuites()
    runner = unittest.TextTestRunner()
    runner.run(suite)
        
if __name__ == '__main__':
    main()
                
