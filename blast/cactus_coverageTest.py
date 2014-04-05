import unittest, os
from sonLib.bioio import getTempFile, popenCatch, system
from textwrap import dedent
from random import sample
from cactus.shared.test import getCactusInputs_blanchette

class TestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        # simple test data -- not an actual alignment, but to test if
        # coverage is correct. no overlap on B, but overlap on A.
        self.simpleFastaPathA = getTempFile()
        open(self.simpleFastaPathA, 'w').write(dedent('''\
        >simpleSeqA1
        ACTAGAGTAGGAGAGAGAGGGGGG
        CATGCATGCATGCATGCATGCATG
        >simpleSeqA2
        AAAAAAAAAAAAAAAACTCGTGAG
        CATGCATGCATGCATGCATGCATG'''))
        self.simpleFastaPathB = getTempFile()
        open(self.simpleFastaPathB, 'w').write(dedent('''\
        >simpleSeqB1
        CATGCATGCATGCATGCATGCATG
        CATGCATGCATGCATGCATGCATG'''))
        self.simpleCigarPath = getTempFile()
        open(self.simpleCigarPath, 'w').write(dedent('''\
        cigar: simpleSeqB1 0 9 + simpleSeqA1 10 0 - 0 M 8 D 1 M 1
        cigar: simpleSeqB1 9 18 + simpleSeqA1 2 6 + 0 M 3 I 5 M 1
        cigar: simpleSeqB1 18 28 + simpleSeqA2 0 10 + 0 M 1 I 2 M 2 D 2 M 5
        cigar: simpleSeqB1 28 30 + simpleSeqA2 6 8 + 0 M 2
        cigar: simpleSeqB1 30 32 + simpleSeqA2 7 9 + 0 M 2
        '''))

    def tearDown(self):
        unittest.TestCase.tearDown(self)
        os.remove(self.simpleFastaPathA)
        os.remove(self.simpleFastaPathB)
        os.remove(self.simpleCigarPath)

    def testSimpleCoverageOnA(self):
        # Genome A
        bed = popenCatch("cactus_coverage %s %s" % (self.simpleFastaPathA,
                                                    self.simpleCigarPath))
        self.assertEqual(bed, dedent('''\
        simpleSeqA1\t0\t1\t\t1
        simpleSeqA1\t2\t5\t\t2
        simpleSeqA1\t5\t6\t\t2
        simpleSeqA1\t6\t10\t\t1
        simpleSeqA2\t0\t1\t\t1
        simpleSeqA2\t1\t3\t\t1
        simpleSeqA2\t5\t6\t\t1
        simpleSeqA2\t6\t7\t\t2
        simpleSeqA2\t7\t8\t\t3
        simpleSeqA2\t8\t9\t\t2
        simpleSeqA2\t9\t10\t\t1
        '''))

    def testSimpleCoverageOnB(self):
        # Genome B
        bed = popenCatch("cactus_coverage %s %s" % (self.simpleFastaPathB,
                                                    self.simpleCigarPath))
        self.assertEqual(bed, dedent('''\
        simpleSeqB1\t0\t8\t\t1
        simpleSeqB1\t8\t9\t\t1
        simpleSeqB1\t9\t12\t\t1
        simpleSeqB1\t17\t18\t\t1
        simpleSeqB1\t18\t19\t\t1
        simpleSeqB1\t21\t23\t\t1
        simpleSeqB1\t23\t28\t\t1
        simpleSeqB1\t28\t30\t\t1
        simpleSeqB1\t30\t32\t\t1
        '''))

    def testInvariants(self):
        (seqs, _) = getCactusInputs_blanchette()
        seqs = sample(seqs, 2)
        cigarPath = getTempFile()
        system("cactus_lastz --format=cigar %s[multiple] %s[multiple] > %s" % \
               (seqs[0], seqs[1], cigarPath))
        bed = popenCatch("cactus_coverage %s %s" % (seqs[0], cigarPath))
        prevChrom = None
        prevStart = None
        prevEnd = None
        # Check that everything is sorted and there are no overlaps
        for line in bed.split("\n"):
            line.strip()
            if line == "":
                continue
            fields = line.split()
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            if chrom == prevChrom:
                self.assertTrue(start > prevStart)
                self.assertTrue(start >= prevEnd)
        os.remove(cigarPath)

if __name__ == '__main__':
    unittest.main()
