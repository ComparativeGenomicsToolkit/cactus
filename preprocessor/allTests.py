#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
import unittest

from cactus.preprocessor.lastzRepeatMasking.cactus_lastzRepeatMaskTest import TestCase as repeatMaskTest
from cactus.preprocessor.cactus_preprocessorTest import TestCase as preprocessorTest

def allSuites():
    allTests = unittest.TestSuite((unittest.makeSuite(repeatMaskTest, 'test'),
                                   unittest.makeSuite(preprocessorTest, 'test')))
    return allTests

def main():
    suite = allSuites()
    runner = unittest.TextTestRunner()
    i = runner.run(suite)
    return len(i.failures) + len(i.errors)
        
if __name__ == '__main__':
    import sys
    sys.exit(main())
