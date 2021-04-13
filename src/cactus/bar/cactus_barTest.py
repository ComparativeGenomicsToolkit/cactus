#!/usr/bin/env python3

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import os
import sys
import random

from sonLib.bioio import TestStatus
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import fastaAlignmentRead
from sonLib.bioio import fastaWrite
from sonLib.bioio import fastaAlignmentWrite
from sonLib.bioio import fastaReadHeaders
from sonLib.bioio import getLogLevelString

from cactus.shared.common import cactus_call

"""Tests cactus_bar. Requires the installation of cactusTools and mafTools.
"""

class TestCase(unittest.TestCase):

    def setUp(self):
        self.testNo = TestStatus.getTestSetup(3, 10, 0, 0)
        self.batchSystem = "parasol"
        unittest.TestCase.setUp(self)

    @TestStatus.shortLength
    def testPosetAlignerAPI(self):
        """Run all the cactus base aligner CuTests, fail if any of them fail.
        """
        cactus_call(parameters=["cactus_barTests", getLogLevelString()])

if __name__ == '__main__':
    unittest.main()
