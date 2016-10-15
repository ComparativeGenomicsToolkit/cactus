#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys

from sonLib.bioio import TestStatus

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_blanchette
from cactus.shared.test import runWorkflow_multipleExamples
from cactus.shared.test import silentOnSuccess

class TestCase(unittest.TestCase):
    @silentOnSuccess
    def testCactus_Random(self):
        runWorkflow_multipleExamples(getCactusInputs_random, 
                                     testNumber=TestStatus.getTestSetup(), 
                                     buildAvgs=True)
    @silentOnSuccess
    def testCactus_Blanchette(self):
        runWorkflow_multipleExamples(getCactusInputs_blanchette, 
                                     testRestrictions=(TestStatus.TEST_SHORT,), inverseTestRestrictions=True, 
                                     buildAvgs=True)
    
if __name__ == '__main__':
    unittest.main()
