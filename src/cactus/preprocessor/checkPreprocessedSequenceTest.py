import os
import random
import shutil
import tempfile
import unittest

from cactus.preprocessor.checkPreprocessedSequence import check_sequence_preserved
from cactus.preprocessor.checkPreprocessedSequence import fasta_digest, hardmask_ok

"""Tests the guard that stops a preprocessor from quietly altering the sequence.

Renaming, soft-masking, re-wrapping and reordering are all things the preprocessors
legitimately do and must be accepted.  A changed, truncated, dropped or added
sequence is data corruption and must be fatal.
"""


class TestCase(unittest.TestCase):
    def setUp(self):
        self.tempDir = tempfile.mkdtemp()
        random.seed(1)
        self.records = [('chr1', self.randomSequence(5000)),
                        ('chr2', self.randomSequence(3000)),
                        ('chr3', self.randomSequence(700))]
        self.inputPath = self.writeFasta('in.fa', self.records)

    def tearDown(self):
        shutil.rmtree(self.tempDir, ignore_errors=True)

    def randomSequence(self, length):
        return ''.join(random.choice('ACGT') for _ in range(length))

    def writeFasta(self, name, records, wrap=60):
        path = os.path.join(self.tempDir, name)
        with open(path, 'w') as fh:
            for header, sequence in records:
                fh.write('>%s\n' % header)
                for i in range(0, len(sequence), wrap):
                    fh.write(sequence[i:i + wrap] + '\n')
        return path

    def check(self, name, records, wrap=60, allow_hardmask=False):
        """Run the guard against a fasta built from records, returning the error text."""
        path = self.writeFasta(name, records, wrap)
        try:
            check_sequence_preserved(self.inputPath, path, event_name='testGenome',
                                     step_name='red', allow_hardmask=allow_hardmask)
        except RuntimeError as e:
            return str(e)
        return None

    # things the preprocessors are allowed to do

    def testIdentical(self):
        self.assertIsNone(self.check('same.fa', self.records))

    def testRenamed(self):
        renamed = [('renamed%d' % i, s) for i, (_, s) in enumerate(self.records)]
        self.assertIsNone(self.check('renamed.fa', renamed))

    def testSoftMasked(self):
        masked = [(n, s[:100].lower() + s[100:]) for n, s in self.records]
        self.assertIsNone(self.check('masked.fa', masked))

    def testRewrapped(self):
        self.assertIsNone(self.check('wrap80.fa', self.records, wrap=80))
        self.assertIsNone(self.check('wrap7.fa', self.records, wrap=7))

    def testReordered(self):
        self.assertIsNone(self.check('reordered.fa', list(reversed(self.records))))

    def testRenamedReorderedAndMasked(self):
        mixed = [('x', self.records[2][1].lower()),
                 ('y', self.records[0][1]),
                 ('z', self.records[1][1])]
        self.assertIsNone(self.check('mixed.fa', mixed))

    # things that are corruption

    def testTruncatedByOneBase(self):
        cut = [(n, s[:-1]) if n == 'chr2' else (n, s) for n, s in self.records]
        self.assertIn('changed the length', self.check('cut1.fa', cut))

    def testTruncatedHeavily(self):
        cut = [(n, s[:100]) if n == 'chr1' else (n, s) for n, s in self.records]
        error = self.check('cut2.fa', cut)
        self.assertIn('changed the length', error)
        self.assertIn('chr1', error)

    def testSequenceDropped(self):
        self.assertIn('changed the number', self.check('dropped.fa', self.records[:2]))

    def testSequenceAdded(self):
        extra = self.records + [('chr4', self.randomSequence(400))]
        self.assertIn('changed the number', self.check('added.fa', extra))

    def testBaseSubstituted(self):
        def substitute(sequence):
            replacement = 'T' if sequence[10] != 'T' else 'A'
            return sequence[:10] + replacement + sequence[11:]
        changed = [(n, substitute(s)) if n == 'chr3' else (n, s) for n, s in self.records]
        self.assertIn('changed the bases', self.check('substituted.fa', changed))

    # hard-masking, which is allowed only where it is configured

    def testHardMaskAllowed(self):
        hardMasked = [(n, s[:50] + 'N' * 200 + s[250:]) for n, s in self.records]
        self.assertIsNone(self.check('hard.fa', hardMasked, allow_hardmask=True))
        self.assertIn('changed the bases', self.check('hard2.fa', hardMasked))

    def testHardMaskDoesNotExcuseOtherChanges(self):
        altered = [(n, s[:50] + 'N' * 200 + 'ACGT' + s[254:]) if n == 'chr1' else (n, s)
                   for n, s in self.records]
        self.assertIn('not hard-masking',
                      self.check('hard3.fa', altered, allow_hardmask=True))

    def testHardMaskOkHelper(self):
        self.assertTrue(hardmask_ok(b'ACGTACGT', b'ACGTACGT'))
        self.assertTrue(hardmask_ok(b'ACGTACGT', b'ACNNACGT'))
        self.assertFalse(hardmask_ok(b'ACGTACGT', b'ACGTACGA'))
        self.assertFalse(hardmask_ok(b'ACGTACGT', b'ACGTACG'))

    # the digest itself

    def testDigestIgnoresCaseAndWrapping(self):
        wrapped = self.writeFasta('w.fa', self.records, wrap=13)
        masked = self.writeFasta('m.fa', [(n, s.lower()) for n, s in self.records])
        expected = [(l, c) for _, l, c in fasta_digest(self.inputPath)]
        self.assertEqual(expected, [(l, c) for _, l, c in fasta_digest(wrapped)])
        self.assertEqual(expected, [(l, c) for _, l, c in fasta_digest(masked)])

    def testDigestLengths(self):
        self.assertEqual([('chr1', 5000), ('chr2', 3000), ('chr3', 700)],
                         [(n, l) for n, l, _ in fasta_digest(self.inputPath)])


if __name__ == '__main__':
    unittest.main()
