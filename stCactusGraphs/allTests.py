#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys
from cactus.shared.test import parseCactusSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import system
from sonLib.bioio import getLogLevelString

class TestCase(unittest.TestCase):
    def testCuTest(self):
        system("stPinchGraphTests %s" % getLogLevelString())
    
def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()