import unittest
from io import StringIO
from textwrap import dedent
from sonLib.bioio import getTempFile
from sonLib.bioio import TestStatus
from cactus.blast.trimSequences import trimSequences
import os

class TestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.faPath = getTempFile()
        open(self.faPath, 'w').write(dedent('''\
        >seq1
        CATGCATGCATGCATGCATGCATG
        CATGCATGCATGCATGCATGCATG
        CATGCATGCATGCATGCATGCATG
        CATGCATGCATGCATGCATGCATG
        >seq2
        ACTGACTGACTGACTGACTGACTG
        ACTGACTGACTGACTGACTGACTG
        ACTGACTGACTGACTGACTGACTG
        ACTGACTGACTGACTGACTGACTG'''))
        self.bedPath = getTempFile()
        open(self.bedPath, 'w').write(dedent('''\
        seq1\t0\t5\t\t1
        seq1\t6\t11\t\t2
        seq1\t15\t16\t\t5'''))

    def tearDown(self):
        os.remove(self.faPath)
        os.remove(self.bedPath)

    @TestStatus.shortLength
    def testSimplestParameters(self):
        # Test w/ no windowing, minimum size, etc to see if bed
        # import/fasta export works
        output = StringIO()
        trimSequences(self.faPath, self.bedPath, output, flanking=0, minSize=0, windowSize=1, threshold=1)
        self.assertTrue(dedent('''\
        >seq1|0
        CATGC''') in output.getvalue())
        self.assertTrue(dedent('''\
        >seq1|6
        TGCAT''') in output.getvalue())
        self.assertTrue(dedent('''\
        >seq1|15
        G''') in output.getvalue())

    @TestStatus.shortLength
    def testComplement(self):
        output = StringIO()
        trimSequences(self.faPath, self.bedPath, output, flanking=0, minSize=0, windowSize=1, threshold=1,
                      complement=True)
        self.assertTrue(dedent('''\
        >seq1|5
        A''') in output.getvalue())
        self.assertTrue(dedent('''\
        >seq1|11''') in output.getvalue())
        self.assertTrue(dedent('''\
        >seq1|16''') in output.getvalue())
        # make sure the sequence that isn't covered at all is included
        self.assertTrue(dedent('''\
        >seq2|0''') in output.getvalue())

    @TestStatus.shortLength
    def testFlanking(self):
        output = StringIO()
        trimSequences(self.faPath, self.bedPath, output, flanking=1, minSize=0, windowSize=1, threshold=1)
        # The two blocks 0-5, 6-11 should be merged together since
        # their flanking sequence intersects. Additionally the
        # flanking sequence shouldn't go past the beginning sequence.
        self.assertTrue(dedent('''\
        >seq1|0
        CATGCATGCATG''') in output.getvalue())
        self.assertTrue(dedent('''\
        >seq1|14
        TGC''') in output.getvalue())

    @TestStatus.shortLength
    def testDepth(self):
        output = StringIO()
        trimSequences(self.faPath, self.bedPath, output, flanking=0, minSize=0, windowSize=1, depth=2)
        self.assertTrue(">seq1|0" not in output.getvalue())
        self.assertTrue(">seq1|6" in output.getvalue())
        self.assertTrue(">seq1|15" in output.getvalue())

    @TestStatus.shortLength
    def testMinSize(self):
        output = StringIO()
        trimSequences(self.faPath, self.bedPath, output, flanking=0, minSize=2, windowSize=1, threshold=1)
        self.assertTrue(">seq1|0" in output.getvalue())
        self.assertTrue(">seq1|6" in output.getvalue())
        self.assertTrue(">seq1|15" not in output.getvalue())

    @TestStatus.shortLength
    def testWithBlankLines(self):
        output = StringIO()
        with open(self.faPath, 'a') as f:
            f.write("\n\n\n")
        trimSequences(self.faPath, self.bedPath, output, flanking=0, minSize=0, windowSize=1, threshold=1)
        self.assertTrue(dedent('''\
        >seq1|0
        CATGC''') in output.getvalue())
        self.assertTrue(dedent('''\
        >seq1|6
        TGCAT''') in output.getvalue())
        self.assertTrue(dedent('''\
        >seq1|15
        G''') in output.getvalue())

if __name__ == "__main__":
    unittest.main()
