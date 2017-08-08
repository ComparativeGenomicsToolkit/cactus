import unittest
from textwrap import dedent
from sonLib.bioio import getTempFile
from cactus.shared.common import cactus_call
from cactus.shared.test import silentOnSuccess
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

    @silentOnSuccess
    def testSimplestParameters(self):
        # Test w/ no windowing, minimum size, etc to see if bed
        # import/fasta export works
        fa = cactus_call(parameters=["cactus_trimSequences.py", "--flanking", "0", "--minSize", "0",
                                     "--windowSize", "1", "--threshold", "1", self.faPath,
                                     self.bedPath],
                         check_output=True)
        self.assertTrue(dedent('''\
        >seq1|0
        CATGC''') in fa)
        self.assertTrue(dedent('''\
        >seq1|6
        TGCAT''') in fa)
        self.assertTrue(dedent('''\
        >seq1|15
        G''') in fa)

    @silentOnSuccess
    def testComplement(self):
        fa = cactus_call(parameters=["cactus_trimSequences.py", "--flanking", "0", "--minSize", "0",
                                     "--windowSize", "1", "--threshold", "1", "--complement",
                                     self.faPath, self.bedPath],
                         check_output=True)
        self.assertTrue(dedent('''\
        >seq1|5
        A''') in fa)
        self.assertTrue(dedent('''\
        >seq1|11''') in fa)
        self.assertTrue(dedent('''\
        >seq1|16''') in fa)
        # make sure the sequence that isn't covered at all is included
        self.assertTrue(dedent('''\
        >seq2|0''') in fa)

    @silentOnSuccess
    def testFlanking(self):
        fa = cactus_call(parameters=["cactus_trimSequences.py", "--flanking", "1", "--minSize", "0",
                                     "--windowSize", "1", "--threshold", "1", self.faPath,
                                     self.bedPath],
                         check_output=True)
        # The two blocks 0-5, 6-11 should be merged together since
        # their flanking sequence intersects. Additionally the
        # flanking sequence shouldn't go past the beginning sequence.
        self.assertTrue(dedent('''\
        >seq1|0
        CATGCATGCATG''') in fa)
        self.assertTrue(dedent('''\
        >seq1|14
        TGC''') in fa)

    @silentOnSuccess
    def testDepth(self):
        fa = cactus_call(parameters=["cactus_trimSequences.py", "--flanking", "0", "--minSize", "0",
                                     "--windowSize", "1", "--depth", "2", self.faPath,
                                     self.bedPath],
                         check_output=True)
        self.assertTrue(">seq1|0" not in fa)
        self.assertTrue(">seq1|6" in fa)
        self.assertTrue(">seq1|15" in fa)

    @silentOnSuccess
    def testMinSize(self):
        fa = cactus_call(parameters=["cactus_trimSequences.py", "--flanking", "0", "--minSize", "2",
                                     "--windowSize", "1", "--threshold", "1", self.faPath,
                                     self.bedPath],
                         check_output=True)
        self.assertTrue(">seq1|0" in fa)
        self.assertTrue(">seq1|6" in fa)
        self.assertTrue(">seq1|15" not in fa)

    @silentOnSuccess
    def testWithBlankLines(self):
        with open(self.faPath, 'a') as f:
            f.write("\n\n\n")
        fa = cactus_call(parameters=["cactus_trimSequences.py", "--flanking", "0", "--minSize", "0",
                                     "--windowSize", "1", "--threshold", "1", self.faPath,
                                     self.bedPath],
                         check_output=True)
        self.assertTrue(dedent('''\
        >seq1|0
        CATGC''') in fa)
        self.assertTrue(dedent('''\
        >seq1|6
        TGCAT''') in fa)
        self.assertTrue(dedent('''\
        >seq1|15
        G''') in fa)

if __name__ == "__main__":
    unittest.main()
