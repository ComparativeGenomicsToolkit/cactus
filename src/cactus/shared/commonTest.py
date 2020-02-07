import os
import shutil
import unittest
from base64 import b64encode

from sonLib.bioio import TestStatus
from sonLib.bioio import getTempFile
from sonLib.bioio import getTempDirectory
from sonLib.bioio import system
from toil.job import Job
from toil.common import Toil
from cactus.shared.common import encodeFlowerNames, decodeFirstFlowerName, \
                                 runCactusSplitFlowersBySecondaryGrouping, \
                                 cactus_call, ChildTreeJob

class TestCase(unittest.TestCase):
    def setUp(self):
        self.testNo = TestStatus.getTestSetup(1, 5, 10, 100)
        self.tempDir = getTempDirectory(os.getcwd())
        self.tempFiles = []
        unittest.TestCase.setUp(self)

    def tearDown(self):
        unittest.TestCase.tearDown(self)
        system("rm -rf %s" % self.tempDir)

    @TestStatus.shortLength
    def testEncodeFlowerNames(self):
        self.assertEqual("3 100 -95 995", encodeFlowerNames([ 100, 5, 1000 ]))
        self.assertEqual("0", encodeFlowerNames([  ]))
        self.assertEqual("1 1", encodeFlowerNames([ 1 ]))

    @TestStatus.shortLength
    def testDecodeFirstFlowerName(self):
        self.assertEqual(None, decodeFirstFlowerName("0 b"))
        self.assertEqual(None, decodeFirstFlowerName("0"))
        self.assertEqual(-1, decodeFirstFlowerName("1 b -1"))
        self.assertEqual(1, decodeFirstFlowerName("2 1 a 1"))
        self.assertEqual(3, decodeFirstFlowerName("2 3 a 1"))
        self.assertEqual(5, decodeFirstFlowerName("2 5 1"))
        self.assertEqual(7, decodeFirstFlowerName("2 b 7 a 1"))
        self.assertEqual(9, decodeFirstFlowerName("4 9 1 1 b 1"))
        self.assertEqual(13, decodeFirstFlowerName("1 b 13"))

    @TestStatus.shortLength
    def testRunCactusSplitFlowersBySecondaryGrouping(self):
        self.assertEqual([(True, "1 -1") ], runCactusSplitFlowersBySecondaryGrouping("1 b -1"))
        self.assertEqual([(False, "1 1"), (False, "1 2")], runCactusSplitFlowersBySecondaryGrouping("2 1 a 1"))
        self.assertEqual([(False, "1 3"), (False, "1 4")], runCactusSplitFlowersBySecondaryGrouping("2 3 a 1"))
        self.assertEqual([(False, "2 5 1")], runCactusSplitFlowersBySecondaryGrouping("2 5 1"))
        self.assertEqual([(True, "1 7"), (False, "1 8")], runCactusSplitFlowersBySecondaryGrouping("2 b 7 a 1"))
        self.assertEqual([(False, "3 9 1 1"), (True, "1 12")], runCactusSplitFlowersBySecondaryGrouping("4 9 1 1 b 1"))
        self.assertEqual([(True, "1 13") ], runCactusSplitFlowersBySecondaryGrouping("1 b 13"))
        self.assertEqual([(False, "3 9 1 1"), (False, "2 8 4"), (True, "3 13 7 8")], runCactusSplitFlowersBySecondaryGrouping("8 9 1 1 a -3 4 b 1 7 8"))

    @TestStatus.shortLength
    def testCactusCall(self):
        inputFile = getTempFile(rootDir=self.tempDir)

        with open("/dev/urandom", "rb") as randText:
            with open(inputFile, 'w') as fh:
                fh.write(b64encode(randText.read(1024)).decode())
        input = "".join(open(inputFile).read().split("\n"))

        #Send input to container's stdin through a file, get output
        #from stdout
        output = "".join(cactus_call(infile=inputFile, check_output=True,
                                     parameters=["docker_test_script"]).split("\n"))
        self.assertEqual(input, output)


        #Send input as string, get output from stdout
        output = "".join(cactus_call(stdin_string=input, check_output=True,
                             parameters=["docker_test_script"]).split("\n"))

        self.assertEqual(input, output)

    @TestStatus.shortLength
    def testCactusCallPipes(self):
        inputFile = getTempFile(rootDir=self.tempDir)
        with open(inputFile, 'w') as f:
            f.write('foobar\n')
        # using 'cat' here rather than infile is intentional; it tests
        # whether the directory is mounted into containers correctly.
        output = cactus_call(parameters=[['cat', inputFile],
                                         ['sed', 's/foo/baz/g'],
                                         ['awk', '{ print "quux" $0 }']],
                             check_output=True)
        self.assertEqual(output, 'quuxbazbar\n')

    @TestStatus.mediumLength
    def testChildTreeJob(self):
        """Check that the ChildTreeJob class runs all children."""
        numChildren = 100
        flagDir = getTempDirectory()

        options = Job.Runner.getDefaultOptions(getTempDirectory())
        shutil.rmtree(options.jobStore)

        with Toil(options) as toil:
            toil.start(CTTestParent(flagDir, numChildren))

        # Check that all jobs ran
        for i in range(numChildren):
            self.assertTrue(os.path.exists(os.path.join(flagDir, str(i))))
        shutil.rmtree(flagDir)

class CTTestParent(ChildTreeJob):
    def __init__(self, flagDir, numChildren):
        self.flagDir = flagDir
        self.numChildren = numChildren
        super(CTTestParent, self).__init__()

    def run(self, fileStore):
        for i in range(self.numChildren):
            self.addChild(CTTestChild(self.flagDir, i))

class CTTestChild(Job):
    def __init__(self, flagDir, index):
        self.flagDir = flagDir
        self.index = index
        super(CTTestChild, self).__init__()

    def run(self, fileStore):
        # Mark that this job has run using a flag file
        path = os.path.join(self.flagDir, str(self.index))
        with open(path, 'w') as f:
            # Empty file
            f.write('')

if __name__ == '__main__':
    unittest.main()
