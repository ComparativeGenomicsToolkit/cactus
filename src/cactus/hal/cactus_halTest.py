#!/usr/bin/env python3

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys

from sonLib.bioio import TestStatus, system, getLogLevelString

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import runWorkflow_multipleExamples
from cactus.shared.test import silentOnSuccess

from cactus.shared.common import cactus_call

class TestCase(unittest.TestCase):
    @silentOnSuccess
    @TestStatus.mediumLength
    def testCactusRecursiveHalGenerator_Random(self):
        runWorkflow_multipleExamples(getCactusInputs_random,
                                     testNumber=TestStatus.getTestSetup(),
                                     buildReference=True, buildHal=True, buildFasta=True)

    @silentOnSuccess
    @TestStatus.mediumLength
    def testCactusRecursiveHalGenerator_Blanchette(self):
        runWorkflow_multipleExamples(getCactusInputs_blanchette,
                                     buildReference=True, buildHal=True, buildFasta=True)

    def testHalGeneratorFunctions(self):
        """Run all the CuTests, fail if any of them fail.
        """
        cactus_call(parameters=["cactus_halGeneratorTests", getLogLevelString()])

if __name__ == '__main__':
    unittest.main()
