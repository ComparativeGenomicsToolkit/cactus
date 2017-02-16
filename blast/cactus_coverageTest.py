import unittest, os, random
from sonLib.bioio import getTempFile, popenCatch, system
from textwrap import dedent
from cactus.shared.test import getCactusInputs_encode

class TestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        # simple test data -- not an actual alignment, but to test if
        # coverage is correct. no overlap on B, but overlap on A.
        self.simpleFastaPathA = getTempFile()
        open(self.simpleFastaPathA, 'w').write(dedent('''\
        >id=0|simpleSeqA1 otherTokens thatDon'tMatter
        ACTAGAGTAGGAGAGAGAGGGGGG
        CATGCATGCATGCATGCATGCATG
        >id=1|simpleSeqA2 otherTokens thatDon'tMatter
        AAAAAAAAAAAAAAAACTCGTGAG
        CATGCATGCATGCATGCATGCATG'''))
        self.simpleFastaPathB = getTempFile()
        open(self.simpleFastaPathB, 'w').write(dedent('''\
        >id=2|simpleSeqB1 otherTokens
        CATGCATGCATGCATGCATGCATG
        CATGCATGCATGCATGCATGCATG'''))
        self.simpleFastaPathC = getTempFile()
        open(self.simpleFastaPathC, 'w').write(dedent('''\
        >id=3|simpleSeqC1 otherTokens thatDon'tMatter
        CATGCATGCATGCATGCATGCATG
        CATGCATGCATGCATGCATGCATG'''))
        self.simpleFastaPathD = getTempFile()
        open(self.simpleFastaPathD, 'w').write(dedent('''\
        >id=4|simpleSeqD otherTokens thatDon'tMatter
        CATGCATGCATGCATGCATGCATG
        CATGCATGCATGCATGCATGCATG'''))
        self.simpleCigarPath = getTempFile()
        open(self.simpleCigarPath, 'w').write(dedent('''\
        cigar: id=2|simpleSeqB1 0 9 + id=0|simpleSeqA1 10 0 - 0 M 8 D 1 M 1
        cigar: id=2|simpleSeqB1 9 18 + id=0|simpleSeqA1 2 6 + 0 M 3 I 5 M 1
        cigar: id=2|simpleSeqB1 18 28 + id=1|simpleSeqA2 0 10 + 0 M 1 I 2 M 2 D 2 M 5
        cigar: id=2|simpleSeqB1 28 30 + id=1|simpleSeqA2 6 8 + 0 M 2
        cigar: id=2|simpleSeqB1 30 32 + id=1|simpleSeqA2 7 9 + 0 M 2
        cigar: id=12|simpleSeqZ1 0 1 + id=0|simpleSeqA1 6 7 + 0 M 1
        cigar: id=3|simpleSeqC1 0 5 + id=4|simpleSeqD 0 5 + 0 M 5
        cigar: id=4|simpleSeqD 5 10 + id=3|simpleSeqC1 5 10 + 0 M 5
        cigar: id=3|simpleSeqC1 10 15 + id=3|simpleSeqC1 15 20 + 0 M 5
        cigar: id=303|simpleSeqNonExistent 0 10 + id=3|simpleSeqC1 0 10 + 0 M 10
        '''))

    def tearDown(self):
        unittest.TestCase.tearDown(self)
        os.remove(self.simpleFastaPathA)
        os.remove(self.simpleFastaPathB)
        os.remove(self.simpleFastaPathC)
        os.remove(self.simpleFastaPathD)
        os.remove(self.simpleCigarPath)

    def testSimpleCoverageOnA(self):
        # Genome A
        bed = popenCatch("cactus_coverage %s %s" % (self.simpleFastaPathA,
                                                    self.simpleCigarPath))
        self.assertEqual(bed, dedent('''\
        id=0|simpleSeqA1\t0\t1\t\t1
        id=0|simpleSeqA1\t2\t7\t\t2
        id=0|simpleSeqA1\t7\t10\t\t1
        id=1|simpleSeqA2\t0\t3\t\t1
        id=1|simpleSeqA2\t5\t6\t\t1
        id=1|simpleSeqA2\t6\t7\t\t2
        id=1|simpleSeqA2\t7\t8\t\t3
        id=1|simpleSeqA2\t8\t9\t\t2
        id=1|simpleSeqA2\t9\t10\t\t1
        '''))

    def testSimpleCoverageOnB(self):
        # Genome B
        bed = popenCatch("cactus_coverage %s %s" % (self.simpleFastaPathB,
                                                    self.simpleCigarPath))
        self.assertEqual(bed, dedent('''\
        id=2|simpleSeqB1\t0\t12\t\t1
        id=2|simpleSeqB1\t17\t19\t\t1
        id=2|simpleSeqB1\t21\t32\t\t1
        '''))

    def testDepthByIDOnA(self):
        # Genome A using depthByID: all depths should be 1 except
        # where 2 different id= prefixes align to the same place:
        # position 6
        bed = popenCatch("cactus_coverage --depthById %s %s" % (
            self.simpleFastaPathA, self.simpleCigarPath))
        self.assertEqual(bed, dedent('''\
        id=0|simpleSeqA1\t0\t1\t\t1
        id=0|simpleSeqA1\t2\t6\t\t1
        id=0|simpleSeqA1\t6\t7\t\t2
        id=0|simpleSeqA1\t7\t10\t\t1
        id=1|simpleSeqA2\t0\t3\t\t1
        id=1|simpleSeqA2\t5\t10\t\t1
        '''))

    def testDepthByIDOnB(self):
        # Genome B using depthByID: should be the same as normal
        # except for 30-31, where it should be 2
        bed = popenCatch("cactus_coverage --depthById %s %s" % (
            self.simpleFastaPathB, self.simpleCigarPath))
        self.assertEqual(bed, dedent('''\
        id=2|simpleSeqB1\t0\t12\t\t1
        id=2|simpleSeqB1\t17\t19\t\t1
        id=2|simpleSeqB1\t21\t32\t\t1
        '''))

    def testFromC(self):
        # Test "--from" filtering by filtering for only alignments
        # from/to D on C.
        bed = popenCatch("cactus_coverage %s %s --from %s" % (self.simpleFastaPathC, self.simpleCigarPath, self.simpleFastaPathD))
        self.assertEqual(bed, dedent('''\
        id=3|simpleSeqC1\t0\t10\t\t1
        '''))

    def testInvariants(self):
        if "SON_TRACE_DATASETS" not in os.environ:
            return
        (seqs, _) = getCactusInputs_encode(random.uniform(0, 2))
        # Chimp encode input has duplicate header names.
        seqs = [i for i in seqs if 'chimp' not in i]
        seqs = random.sample(seqs, 2)
        cigarPath = getTempFile()
        system("cPecanLastz --format=cigar %s[multiple] %s[multiple] > %s" % \
               (seqs[0], seqs[1], cigarPath))
        bed = popenCatch("cactus_coverage %s %s" % (seqs[1], cigarPath))
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
            self.assertTrue(end - start >= 1)
            if chrom == prevChrom:
                self.assertTrue(start > prevStart)
                self.assertTrue(start >= prevEnd)
        os.remove(cigarPath)

if __name__ == '__main__':
    unittest.main()
