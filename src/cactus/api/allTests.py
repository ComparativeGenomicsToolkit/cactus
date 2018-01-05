#!/usr/bin/env python2.7

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys
import os
from sonLib.bioio import TestStatus
from sonLib.bioio import system
from sonLib.bioio import getLogLevelString

from cactus.shared.common import cactus_call

class TestCase(unittest.TestCase):

    def setUp(self):
        self.testNo = TestStatus.getTestSetup()
        self.tempFiles = []
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        for tempFile in self.tempFiles:
            os.remove(tempFile)
        unittest.TestCase.tearDown(self)
        
    def testAPI(self):
        """Run all the cactusAPI CuTests, fail if any of them fail.
        """
        cactus_call(parameters=["cactusAPITests", getLogLevelString()])
        
if __name__ == '__main__':
    unittest.main()
