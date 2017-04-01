import unittest
import os

from sonLib.bioio import TestStatus
from sonLib.bioio import getTempFile
from sonLib.bioio import getTempDirectory
from sonLib.bioio import system
from cactus.shared.common import encodeFlowerNames, decodeFirstFlowerName, \
                                 runCactusSplitFlowersBySecondaryGrouping, \
                                 cactus_call


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
        self.assertEquals([(True, "1 -1") ], runCactusSplitFlowersBySecondaryGrouping("1 b -1"))
        self.assertEquals([(False, "1 1"), (False, "1 2")], runCactusSplitFlowersBySecondaryGrouping("2 1 a 1"))
        self.assertEquals([(False, "1 3"), (False, "1 4")], runCactusSplitFlowersBySecondaryGrouping("2 3 a 1"))
        self.assertEquals([(False, "2 5 1")], runCactusSplitFlowersBySecondaryGrouping("2 5 1"))
        self.assertEquals([(True, "1 7"), (False, "1 8")], runCactusSplitFlowersBySecondaryGrouping("2 b 7 a 1"))
        self.assertEquals([(False, "3 9 1 1"), (True, "1 12")], runCactusSplitFlowersBySecondaryGrouping("4 9 1 1 b 1"))
        self.assertEquals([(True, "1 13") ], runCactusSplitFlowersBySecondaryGrouping("1 b 13"))
        self.assertEquals([(False, "3 9 1 1"), (False, "2 8 4"), (True, "3 13 7 8")], runCactusSplitFlowersBySecondaryGrouping("8 9 1 1 a -3 4 b 1 7 8"))

    def testCactusCall(self):
        inputFile = getTempFile(rootDir=self.tempDir)
        
        with open("/dev/urandom") as randText:
            with open(inputFile, 'w') as fh:
                fh.write(randText.read(1024))
        input = "".join(open(inputFile).read().split("\n"))
        
        #Send input to container's stdin through a file, get output
        #from stdout
        output = "".join(cactus_call(infile=inputFile, check_output=True,
                                     parameters=["docker_test_script"]).split("\n"))
        self.assertEquals(input, output)


        #Send input as string, get output from stdout
        output = "".join(cactus_call(stdin_string=input, check_output=True,
                             parameters=["docker_test_script"]).split("\n"))

        self.assertEquals(input, output)
      
if __name__ == '__main__':
    unittest.main()
