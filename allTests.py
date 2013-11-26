#!/usr/bin/env python

#Copyright (C) 2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest

from cactus.setup.cactus_setupTest import TestCase as setupTest
from cactus.blast.cactus_blastTest import TestCase as blastTest
from cactus.blast.cactus_realignTest import TestCase as realignTest
from cactus.pipeline.cactus_workflowTest import TestCase as workflowTest
from cactus.pipeline.cactus_evolverTest import TestCase as evolverTest
from cactus.bar.cactus_barTest import TestCase as barTest
from cactus.phylogeny.cactus_phylogenyTest import TestCase as phylogenyTest
from cactus.faces.cactus_fillAdjacenciesTest import TestCase as adjacenciesTest
from cactus.reference.cactus_referenceTest import TestCase as referenceTest
from cactus.hal.cactus_halTest import TestCase as halTest
from cactus.api.allTests import TestCase as aPITest
from cactus.caf.allTests import TestCase as cafTest
from cactus.normalisation.cactus_normalisationTest import TestCase as normalisationTest
from cactus.progressive.allTests import allSuites as progressiveSuite
from cactus.shared.commonTest import TestCase as commonTest
from cactus.preprocessor.cactus_lastzRepeatMaskTest import TestCase as repeatMaskTest

from cactus.shared.test import parseCactusSuiteTestOptions

def allSuites(): 
    allTests = unittest.TestSuite((unittest.makeSuite(cafTest, 'test'),
                                   unittest.makeSuite(setupTest, 'test'),
                                   unittest.makeSuite(blastTest, 'test'),
                                   unittest.makeSuite(repeatMaskTest, 'test'),
                                   unittest.makeSuite(workflowTest, 'test'),
                                   unittest.makeSuite(evolverTest, 'test'),
                                   unittest.makeSuite(barTest, 'test'),
                                   unittest.makeSuite(realignTest, 'test'),
                                   unittest.makeSuite(phylogenyTest, 'test'),
                                   unittest.makeSuite(adjacenciesTest, 'test'),
                                   unittest.makeSuite(referenceTest, 'test'),
                                   unittest.makeSuite(aPITest, 'test'),
                                   unittest.makeSuite(normalisationTest, 'test'),  
                                   unittest.makeSuite(halTest, 'test'),    
                                   unittest.makeSuite(commonTest, 'test'),                           
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
                
