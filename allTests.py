#!/usr/bin/env python

#Copyright (C) 2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import os
from multiprocessing import Process
import time

from cactus.setup.cactus_setupTest import TestCase as setupTest
from cactus.blast.cactus_blastTest import TestCase as blastTest
from cactus.blast.cactus_coverageTest import TestCase as coverageTest
from cactus.blast.cactus_trimSequencesTest import TestCase as trimSequencesTest
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

def allSuites(): 
    allTests = unittest.TestSuite()
    allTests.addTests([unittest.makeSuite(i) for i in
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
                        commonTest]] + 
                        [progressiveSuite()])
    if "SON_TRACE_DATASETS" in os.environ:
        allTests.addTests([unittest.makeSuite(blastTest), preprocessorTest(), unittest.makeSuite(lastzRepeatMaskTest), unittest.makeSuite(realignTest)])
    return allTests

def main():
    suite = allSuites()
    runner = unittest.TextTestRunner(verbosity=2)
    i = runner.run(suite)
    return len(i.failures) + len(i.errors)

if __name__ == '__main__':
    import sys
    sys.exit(main())
                
