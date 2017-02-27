import unittest
import os
import sys

from sonLib.bioio import TestStatus
from sonLib.bioio import getTempFile
from cactus.shared.common import *

from toil.job import Job


class TestCase(unittest.TestCase):
    def setUp(self):
        self.testNo = TestStatus.getTestSetup(1, 5, 10, 100)
        self.tempDir = getTempDirectory(os.getcwd())
        self.tempFiles = []
        unittest.TestCase.setUp(self)
        
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        system("rm -rf %s" % self.tempDir)
        
    def testEncodeFlowerNames(self):
        self.assertEquals("3 100 -95 995", encodeFlowerNames([ 100, 5, 1000 ]))
        self.assertEquals("0", encodeFlowerNames([  ]))
        self.assertEquals("1 1", encodeFlowerNames([ 1 ]))
    
    def testDecodeFirstFlowerName(self):
        self.assertEquals(None, decodeFirstFlowerName("0 b"))
        self.assertEquals(None, decodeFirstFlowerName("0"))
        self.assertEquals(-1, decodeFirstFlowerName("1 b -1"))
        self.assertEquals(1, decodeFirstFlowerName("2 1 a 1"))
        self.assertEquals(3, decodeFirstFlowerName("2 3 a 1"))
        self.assertEquals(5, decodeFirstFlowerName("2 5 1"))
        self.assertEquals(7, decodeFirstFlowerName("2 b 7 a 1"))
        self.assertEquals(9, decodeFirstFlowerName("4 9 1 1 b 1"))
        self.assertEquals(13, decodeFirstFlowerName("1 b 13"))
    
    def testRunCactusSplitFlowersBySecondaryGrouping(self):
        options = Job.Runner.getDefaultOptions(os.path.join(self.tempDir, "tmpToil"))
        Job.Runner.startToil(Job.wrapJobFn(_testCactusCallFn), options)
    
    @unittest.skip("")
    def testCactusCall(self):
        options = Job.Runner.getDefaultOptions(os.path.join(self.tempDir, "tmpToil"))
        Job.Runner.startToil(Job.wrapJobFn(_testCactusCallFn), options)



def _testRunCactusSplitFlowersBySecondaryGroupingFn(job):
    assert [(True, "1 -1") ] == runCactusSplitFlowersBySecondaryGrouping(job, "1 b -1")
    assert [(False, "1 1"), (False, "1 2")] == runCactusSplitFlowersBySecondaryGrouping(job, "2 1 a 1")
    assert [(False, "1 3"), (False, "1 4")] == runCactusSplitFlowersBySecondaryGrouping(job, "2 3 a 1")
    assert [(False, "2 5 1")] == runCactusSplitFlowersBySecondaryGrouping(job, "2 5 1")
    assert [(True, "1 7"), (False, "1 8")] == runCactusSplitFlowersBySecondaryGrouping(job, "2 b 7 a 1")
    assert [(False, "3 9 1 1"), (True, "1 12")] == runCactusSplitFlowersBySecondaryGrouping(job, "4 9 1 1 b 1")
    assert [(True, "1 13") ] == runCactusSplitFlowersBySecondaryGrouping(job, "1 b 13")
    assert [(False, "3 9 1 1"), (False, "2 8 4"), (True, "3 13 7 8")] == runCactusSplitFlowersBySecondaryGrouping(job, "8 9 1 1 a -3 4 b 1 7 8")

def _testCactusCallFn(job):
    from cactus.shared.common import cactus_call
    inputFile = job.fileStore.getLocalTempFile()

    with open("/dev/urandom") as randText:
        with open(inputFile, 'w') as fh:
            fh.write(randText.read(1024))
    inputString = "".join(open(inputFile).read().split("\n"))

    #Send input to container's stdin through a file, get output
    #from stdout
    output = "".join(cactus_call(job, tool="ubuntu", infile=inputFile, check_output=True,
                                 parameters=["cat"]).split("\n"))
    assert inputString == output


    #Send input as string, get output from stdout
    output = "".join(cactus_call(job, tool="ubuntu", stdin_string=inputString, check_output=True,
                         parameters=["cat"]).split("\n"))

    assert inputString == output
      
if __name__ == '__main__':
    unittest.main()
