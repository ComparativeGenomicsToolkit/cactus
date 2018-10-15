#!/usr/bin/env python

#Copyright (C) 2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import os

from cactus.setup.cactus_setupTest import TestCase as setupTest
from cactus.blast.blastTest import TestCase as blastTest
from cactus.blast.cactus_coverageTest import TestCase as coverageTest
from cactus.blast.trimSequencesTest import TestCase as trimSequencesTest
from cactus.blast.mappingQualityRescoringAndFilteringTest import TestCase as mappingQualityTest
from cactus.pipeline.cactus_workflowTest import TestCase as workflowTest
from cactus.pipeline.cactus_evolverTest import TestCase as evolverTest
from cactus.bar.cactus_barTest import TestCase as barTest
from cactus.phylogeny.cactus_phylogenyTest import TestCase as phylogenyTest
from cactus.faces.cactus_fillAdjacenciesTest import TestCase as adjacenciesTest
from cactus.reference.cactus_referenceTest import TestCase as referenceTest
from cactus.hal.cactus_halTest import TestCase as halTest
from cactus.api.allTests import TestCase as apiTest
from cactus.caf.allTests import TestCase as cafTest
from cactus.normalisation.cactus_normalisationTest import TestCase as normalisationTest
from cactus.progressive.allTests import allSuites as progressiveSuite
from cactus.shared.commonTest import TestCase as commonTest
from cactus.shared.experimentWrapperTest import TestCase as experimentWrapperTest
from cactus.faces.cactus_fillAdjacenciesTest import TestCase as fillAdjacenciesTest
from cactus.preprocessor.allTests import allSuites as preprocessorTest
from cactus.preprocessor.lastzRepeatMasking.cactus_lastzRepeatMaskTest import TestCase as lastzRepeatMaskTest
from cactus.blast.cactus_realignTest import TestCase as realignTest

def getSuite():
    SLOW_BLAST_SUITE = [unittest.makeSuite(blastTest), unittest.makeSuite(mappingQualityTest),
                        preprocessorTest(), unittest.makeSuite(lastzRepeatMaskTest),
                        unittest.makeSuite(realignTest)]
    NORMAL_SUITE = [unittest.makeSuite(i) for i in
                    [setupTest,
                     cafTest,
                     workflowTest,
                     evolverTest,
                     barTest,
                     phylogenyTest,
                     adjacenciesTest,
                     referenceTest,
                     apiTest,
                     normalisationTest,
                     halTest,
                     coverageTest,
                     trimSequencesTest,
                     experimentWrapperTest,
                     fillAdjacenciesTest,
                     commonTest]] + [progressiveSuite()]

    combinedTests = unittest.TestSuite()
    # CI needs to be able to subset tests because they can take over
    # 50 minutes in debug mode. The CACTUS_TEST_CHOICE env variable
    # determines what suites to run.
    suiteChoice = os.environ.get("CACTUS_TEST_CHOICE")
    if suiteChoice is None:
        # Default: all tests. Includes the SON_TRACE_DATASETS reliant
        # tests if there is test data available.
        print "Running all tests."
        combinedTests.addTests(NORMAL_SUITE)
        if "SON_TRACE_DATASETS" in os.environ:
            combinedTests.addTests(SLOW_BLAST_SUITE)
    elif suiteChoice == "normal":
        print "Running only non-blast tests."
        combinedTests.addTests(NORMAL_SUITE)
    elif suiteChoice == "blast":
        print "Running only blast tests."
        combinedTests.addTests(SLOW_BLAST_SUITE)
    else:
        raise RuntimeError("Unknown suite choice: " + suiteChoice)
    return combinedTests

def main():
    suite = getSuite()
    runner = unittest.TextTestRunner(verbosity=2)
    i = runner.run(suite)
    return len(i.failures) + len(i.errors)

if __name__ == '__main__':
    import sys
    sys.exit(main())
                
