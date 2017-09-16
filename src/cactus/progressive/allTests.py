#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
import unittest

from cactus.progressive.multiCactusTreeTest import TestCase as multiCactusTreeTest
from cactus.progressive.outgroupTest import TestCase as outgroupTest
from cactus.progressive.scheduleTest import TestCase as scheduleTest
from cactus.progressive.cactus_progressiveTest import TestCase as cactus_progressiveTest

def allSuites(): 
    allTests = unittest.TestSuite((unittest.makeSuite(multiCactusTreeTest, 'test'),
                                   unittest.makeSuite(outgroupTest, 'test'),
                                   unittest.makeSuite(scheduleTest, 'test'),
                                   unittest.makeSuite(cactus_progressiveTest, 'test')))
    return allTests
        
def main():
    suite = allSuites()
    runner = unittest.TextTestRunner()
    i = runner.run(suite)
    return len(i.failures) + len(i.errors)
        
if __name__ == '__main__':
    import sys
    sys.exit(main())
