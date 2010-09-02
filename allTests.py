import unittest

from cactus.setup.cactus_setupTest import TestCase as setupTest
from cactus.blastAlignment.cactus_alignerTest import TestCase as alignerTest
from cactus.blastAlignment.cactus_batchTest import TestCase as batchTest
from cactus.core.cactus_coreTest import TestCase as coreTest
from cactus.pipeline.cactus_workflowTest import TestCase as workflowTest
from cactus.baseAlignment.cactus_baseAlignerTest import TestCase as baseAlignerTest
from cactus.phylogeny.cactus_phylogenyTest import TestCase as phylogenyTest
from cactus.faces.cactus_fillAdjacenciesTest import TestCase as adjacenciesTest
from cactus.reference.cactus_referenceTest import TestCase as referenceTest
from cactus.api.allTests import TestCase as aPITest
from cactus.normalisation.cactus_normalisationTest import TestCase as normalisationTest
import cactus.utilities.allTests 
 
from cactus.shared.test import parseCactusSuiteTestOptions

def allSuites(): 
    allTests = unittest.TestSuite((unittest.makeSuite(setupTest, 'test'),
                                   unittest.makeSuite(alignerTest, 'test'),
                                   unittest.makeSuite(batchTest, 'test'),
                                   unittest.makeSuite(coreTest, 'test'),
                                   unittest.makeSuite(workflowTest, 'test'),
                                   unittest.makeSuite(baseAlignerTest, 'test'),
                                   unittest.makeSuite(phylogenyTest, 'test'),
                                   unittest.makeSuite(adjacenciesTest, 'test'),
                                   unittest.makeSuite(referenceTest, 'test'),
                                   unittest.makeSuite(aPITest, 'test'),
                                   unittest.makeSuite(normalisationTest, 'test'),
                                   cactus.utilities.allTests.allSuites()))
    return allTests
        
def main():
    parseCactusSuiteTestOptions()
    
    suite = allSuites()
    runner = unittest.TextTestRunner()
    runner.run(suite)
        
if __name__ == '__main__':
    main()
                
