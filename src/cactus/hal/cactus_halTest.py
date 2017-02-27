#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys
import os
import shutil

from sonLib.bioio import TestStatus, system, getLogLevelString, getTempDirectory

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import runWorkflow_multipleExamples
from cactus.shared.test import silentOnSuccess

from cactus.shared.common import cactus_call

from toil.job import Job

class TestCase(unittest.TestCase):

    def setUp(self):
        self.tempDir = getTempDirectory(os.getcwd())
        unittest.TestCase.setUp(self)
        
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        shutil.rmtree(self.tempDir)
        
    @silentOnSuccess
    @unittest.skip("")
    def testCactusRecursiveHalGenerator_Random(self):
        runWorkflow_multipleExamples(getCactusInputs_random,
                                     testNumber=TestStatus.getTestSetup(),
                                     buildReference=True, buildHal=True, buildFasta=True)

    @silentOnSuccess
    @unittest.skip("")
    def testCactusRecursiveHalGenerator_Blanchette(self):
        runWorkflow_multipleExamples(getCactusInputs_blanchette,
                                     testRestrictions=(TestStatus.TEST_SHORT,), inverseTestRestrictions=True,
                                     buildReference=True, buildHal=True, buildFasta=True)
    @silentOnSuccess
    def testHalGeneratorFunctions(self):
        """Run all the CuTests, fail if any of them fail.
        """
        options = Job.Runner.getDefaultOptions(os.path.join(self.tempDir, "tmpToil"))
        Job.Runner.startToil(Job.wrapJobFn(_testHalGeneratorFunctionsFn), options)


def _testHalGeneratorFunctionsFn(job):
    cactus_call(job, parameters=["cactus_halGeneratorTests", getLogLevelString()])    

if __name__ == '__main__':
    unittest.main()
