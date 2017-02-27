#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys
import os
import shutil
from sonLib.bioio import TestStatus
from sonLib.bioio import system
from sonLib.bioio import getLogLevelString
from sonLib.bioio import getTempDirectory

from cactus.shared.common import cactus_call

from toil.job import Job

class TestCase(unittest.TestCase):

    def setUp(self):
        self.testNo = TestStatus.getTestSetup()
        self.tempFiles = []
        self.tempDir = getTempDirectory(os.getcwd())
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        for tempFile in self.tempFiles:
            os.remove(tempFile)
        shutil.rmtree(self.tempDir)
        unittest.TestCase.tearDown(self)
        
    def testAPI(self):
        """Run all the cactusAPI CuTests, fail if any of them fail.
        """
        options = Job.Runner.getDefaultOptions(os.path.join(self.tempDir, "tmpToil"))
        Job.Runner.startToil(Job.wrapJobFn(_testAPIFn), options)

def _testAPIFn(job):
    cactus_call(job, parameters=["cactusAPITests", getLogLevelString()])
        
if __name__ == '__main__':
    unittest.main()
