import gzip
import os
import random
import shutil
import subprocess
import tempfile
import unittest

from cactus.progressive.multiCactusTree import MultiCactusTree
from sonLib.nxnewick import NXNewick

from cactus.hal import cactus_validate
from cactus.preprocessor.checkPreprocessedSequence import fasta_digest
from cactus.hal.cactus_validate import expected_records
from cactus.hal.cactus_validate import hal_sequence_name
from cactus.hal.cactus_validate import looks_like_pangenome
from cactus.hal.cactus_validate import resolve_subpath
from cactus.hal.cactus_validate import validate_genome
from cactus.hal.cactus_validate import fail_hal_validation
from cactus.hal.cactus_validate import validate_hal_export

"""Tests the check that a HAL still holds the sequence that was aligned.

The HAL side is stubbed out: reading a real one needs a real alignment, and what
is worth testing here is everything around it -- replaying the sanitizer's
renaming, folding the ambiguity codes, resolving a clipped contig back to the
range of the input it came from, and telling corruption apart from the
differences a HAL is entitled to have.
"""

# cactus_sanitizeFastaHeaders is what hal_sequence_name() has to agree with, so
# where it is on the PATH the emulation is pinned against it rather than against
# a table of expected answers
SANITIZER = shutil.which('cactus_sanitizeFastaHeaders')


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

    def writeFasta(self, name, records, wrap=60, compress=False):
        path = os.path.join(self.tempDir, name)
        opener = gzip.open if compress else open
        with opener(path, 'wt') as fh:
            for header, sequence in records:
                fh.write('>%s\n' % header)
                for i in range(0, len(sequence), wrap):
                    fh.write(sequence[i:i + wrap] + '\n')
        return path

    def stubGenomes(self, genomes):
        """Make the HAL look like it holds exactly these genomes."""
        original = cactus_validate.hal_genomes
        cactus_validate.hal_genomes = lambda hal_path: list(genomes)
        self.addCleanup(setattr, cactus_validate, 'hal_genomes', original)

    def stubHal(self, records):
        """Make validate_genome() see exactly these (name, sequence) pairs in the HAL."""
        path = self.writeFasta('hal.fa', records)
        original = cactus_validate.hal_records

        def stub(hal_path, genome, work_dir, lengths_only=False, job_memory=None):
            digest = fasta_digest(path)
            return [(n, l, 0) for n, l, _ in digest] if lengths_only else digest

        cactus_validate.hal_records = stub
        self.addCleanup(setattr, cactus_validate, 'hal_records', original)

    def check(self, halRecords, inputPath=None, **kwargs):
        """Validate a stubbed HAL against an input fasta, returning the problems."""
        self.stubHal(halRecords)
        return validate_genome('unused.hal', 'testGenome', inputPath or self.inputPath,
                               self.tempDir, **kwargs)

    # replaying the sanitizer's renaming

    def testNamePlain(self):
        self.assertEqual('chr1', hal_sequence_name('chr1 with a description'))
        self.assertEqual('chr1', hal_sequence_name('id=testGenome|chr1'))
        self.assertEqual('a|pipe|name', hal_sequence_name('a|pipe|name'))
        # only pangenome mode touches '#' and ':'
        self.assertEqual('HG002#1#chr1', hal_sequence_name('HG002#1#chr1'))
        self.assertEqual('chr1:10-15', hal_sequence_name('chr1:10-15'))

    def testNamePangenome(self):
        self.assertEqual('chr1', hal_sequence_name('HG002#1#chr1', pangenome=True))
        self.assertEqual('chr1_sub_9_15', hal_sequence_name('HG002#1#chr1:10-15', pangenome=True))
        # samtools treats a 0 start like a 1, and so does the sanitizer
        self.assertEqual('chr1_sub_0_15', hal_sequence_name('chr1:0-15', pangenome=True))
        self.assertEqual('weird_name_here', hal_sequence_name('weird:name:here', pangenome=True))
        self.assertEqual('chr1', hal_sequence_name('id=HG002.1|chr1', pangenome=True))

    @unittest.skipIf(SANITIZER is None, 'cactus_sanitizeFastaHeaders not on PATH')
    def testNameMatchesTheSanitizer(self):
        headers = ['chr1 with a description', 'HG002#1#chr2:11-20', 'HG002#1#chr3',
                   'id=OTHER|chr4', 'weird:name:here', 'a|pipe|name', 'chr5:0-4']
        path = self.writeFasta('headers.fa', [(h, 'ACGT') for h in headers])
        for pangenome in (False, True):
            command = [SANITIZER, path, 'MYGENOME'] + (['-p'] if pangenome else [])
            output = subprocess.check_output(command, text=True)
            sanitized = []
            for line in output.splitlines():
                if line.startswith('>'):
                    # cactus_consolidated strips the prefix the sanitizer added
                    name = line[1:]
                    if name.startswith('id=') and '|' in name:
                        name = name[name.find('|') + 1:]
                    sanitized.append(name)
            self.assertEqual(sanitized,
                             [hal_sequence_name(h, pangenome=pangenome) for h in headers])

    # digesting the input the way the HAL will hold it

    def testEmptySequencesAreDropped(self):
        path = self.writeFasta('empty.fa', self.records + [('chr4', '')])
        records, _ = expected_records(path)
        self.assertEqual(['chr1', 'chr2', 'chr3'], [n for n, _, _ in records])

    def testAmbiguityCodesFoldToN(self):
        folded = self.writeFasta('iupac.fa', [('chr1', 'ACGTRYKM')])
        ns = self.writeFasta('ns.fa', [('chr1', 'ACGTNNNN')])
        self.assertEqual(expected_records(folded)[0], expected_records(ns)[0])

    def testGzippedInput(self):
        gzipped = self.writeFasta('in.fa.gz', self.records, compress=True)
        self.assertEqual(expected_records(self.inputPath)[0], expected_records(gzipped)[0])

    def testDirectoryOfFastas(self):
        directory = os.path.join(self.tempDir, 'parts')
        os.makedirs(directory)
        for name, sequence in self.records:
            with open(os.path.join(directory, name + '.fa'), 'w') as fh:
                fh.write('>%s\n%s\n' % (name, sequence))
        self.assertEqual(sorted(expected_records(self.inputPath)[0]),
                         sorted(expected_records(directory)[0]))

    # resolving a clipped contig back to the input it came from

    def testSubpathResolution(self):
        lengths = {'chr1': 100}
        self.assertEqual(('chr1', 10, 20), resolve_subpath('chr1_sub_10_20', lengths))
        # clipping can nest, and the inner range is relative to the outer piece
        self.assertEqual(('chr1', 12, 15), resolve_subpath('chr1_sub_10_20_sub_2_5', lengths))
        self.assertIsNone(resolve_subpath('chr1_sub_10_200', lengths))
        self.assertIsNone(resolve_subpath('chr1_sub_20_10', lengths))
        self.assertIsNone(resolve_subpath('chr9_sub_10_20', lengths))
        self.assertIsNone(resolve_subpath('chr1_sub_ten_twenty', lengths))

    # what a HAL is allowed to differ by

    def testIdentical(self):
        self.assertEqual([], self.check(self.records))

    def testSoftMasked(self):
        masked = [(n, s[:100].lower() + s[100:]) for n, s in self.records]
        self.assertEqual([], self.check(masked))

    def testReorderedAndRewrapped(self):
        self.assertEqual([], self.check(list(reversed(self.records))))

    def testPrefixedInput(self):
        prefixed = self.writeFasta('prefixed.fa',
                                   [('id=testGenome|' + n, s) for n, s in self.records])
        self.assertEqual([], self.check(self.records, inputPath=prefixed))

    def testClippedSequence(self):
        name, sequence = self.records[0]
        clipped = [('%s_sub_100_600' % name, sequence[100:600])] + self.records[1:]
        self.assertEqual([], self.check(clipped))

    def testClippedSequenceWithTheWrongOffset(self):
        name, sequence = self.records[0]
        clipped = [('%s_sub_100_600' % name, sequence[101:601])] + self.records[1:]
        problems = self.check(clipped)
        self.assertEqual(1, len(problems))
        self.assertIn('changed the bases', problems[0])

    # what it is not

    def testTruncated(self):
        cut = [(n, s[:-1]) if n == 'chr2' else (n, s) for n, s in self.records]
        problems = self.check(cut)
        self.assertEqual(1, len(problems))
        self.assertIn('changed the length', problems[0])
        self.assertIn('chr2', problems[0])

    def testBaseSubstituted(self):
        def substitute(sequence):
            return sequence[:10] + ('T' if sequence[10] != 'T' else 'A') + sequence[11:]
        changed = [(n, substitute(s)) if n == 'chr3' else (n, s) for n, s in self.records]
        problems = self.check(changed)
        self.assertEqual(1, len(problems))
        self.assertIn('changed the bases', problems[0])

    def testSequenceMissingFromHal(self):
        problems = self.check(self.records[:2])
        self.assertEqual(1, len(problems))
        self.assertIn('is missing 1 of the 3 input sequence(s)', problems[0])
        # --allowMissing keeps the finding but marks it as something to report
        # rather than something to fail on
        relaxed = self.check(self.records[:2], allow_missing=True)
        self.assertEqual(1, len(relaxed))
        self.assertTrue(relaxed[0].startswith('WARNING: '))

    def testSequenceNotInTheInput(self):
        extra = self.records + [('chr4', self.randomSequence(400))]
        problems = self.check(extra)
        self.assertEqual(1, len(problems))
        self.assertIn('not in the input fasta', problems[0])

    def testDuplicateInputNames(self):
        # two headers that sanitise to the same name; the sanitizer refuses these
        path = self.writeFasta('dup.fa', [('chr1 one', 'ACGT'), ('chr1 two', 'ACGT')])
        problems = self.check(self.records, inputPath=path)
        self.assertEqual(1, len(problems))
        self.assertIn('duplicated sequence name', problems[0])

    # working out which naming convention a HAL used

    def testPangenomeShapeDetection(self):
        def tree(newick):
            return MultiCactusTree(NXNewick().parseString(newick, addImpliedRoots=False))
        # a per-chromosome pangenome alignment: star, minigraph leaf still there
        self.assertTrue(looks_like_pangenome(
            tree('(_MINIGRAPH_:1,HG002.1:1,HG002.2:1)Anc0;'), '_MINIGRAPH_'))
        # the merged one: graphmap-join strips the minigraph leaf out
        self.assertTrue(looks_like_pangenome(
            tree('(HG002.1:1,HG002.2:1,CHM13:1)Anc0;'), '_MINIGRAPH_'))
        # a progressive alignment has internal nodes
        self.assertFalse(looks_like_pangenome(
            tree('((a:1,b:1)Anc1:1,(c:1,d:1)Anc2:1)Anc0;'), '_MINIGRAPH_'))
        # ... and every two-genome subproblem cactus-prepare emits is a star
        self.assertFalse(looks_like_pangenome(
            tree('(simMouse_chr6:1,simRat_chr6:1)mr;'), '_MINIGRAPH_'))

    def testConventionFallsBackToWhatTheHalSays(self):
        """A star-shaped progressive alignment gets guessed wrong; the names settle it."""
        pangenome_input = self.writeFasta('pg.fa', [('HG002#1#' + n, s) for n, s in self.records])
        # guessed progressive, but only the pangenome convention reproduces the names
        self.assertEqual([], self.check(self.records, inputPath=pangenome_input,
                                        pangenome=False))
        # guessed pangenome, but these names must be left alone
        ranged = [('chr1:1-100', self.records[0][1])]
        ranged_input = self.writeFasta('ranged.fa', ranged)
        self.assertEqual([], self.check(ranged, inputPath=ranged_input, pangenome=True))

    def testForcedConventionIsNotSecondGuessed(self):
        ranged = [('chr1:1-100', self.records[0][1])]
        ranged_input = self.writeFasta('ranged2.fa', ranged)
        problems = self.check(ranged, inputPath=ranged_input, pangenome=True,
                              force_convention=True)
        self.assertEqual(2, len(problems))
        self.assertIn('not in the input fasta', problems[0])

    # the hook export_hal calls

    def testExportHookAcceptsAndRejects(self):
        self.stubHal(self.records)
        self.stubGenomes(['Anc0', 'testGenome'])
        tree = _FakeTree(['Anc0', 'testGenome'])
        seq_id_map = {'testGenome': self.inputPath}

        job = _FakeJob()
        self.assertEqual([], validate_hal_export(job, 'unused.hal', self.tempDir, tree,
                                                 'Anc0', seq_id_map))
        self.assertTrue(any('Validated' in line for line in job.log))

        self.stubHal([(n, s[:-1]) if n == 'chr2' else (n, s) for n, s in self.records])
        job = _FakeJob()
        problems = validate_hal_export(job, 'unused.hal', self.tempDir, tree, 'Anc0',
                                       seq_id_map)
        self.assertEqual(1, len(problems))
        self.assertIn('changed the length', problems[0])
        # nothing claims to have validated an alignment that did not check out
        self.assertFalse(any('Validated' in line for line in job.log))
        # and the run is failed from a job of its own, so that retrying the
        # failure does not mean rebuilding the alignment
        with self.assertRaises(RuntimeError) as caught:
            fail_hal_validation(job, 'unused.hal', problems)
        self.assertIn('changed the length', str(caught.exception))

    def testExportHookNoticesAMissingGenome(self):
        self.stubHal(self.records)
        self.stubGenomes(['Anc0'])
        problems = validate_hal_export(_FakeJob(), 'unused.hal', self.tempDir,
                                       _FakeTree(['Anc0', 'testGenome']), 'Anc0',
                                       {'testGenome': self.inputPath})
        self.assertEqual(1, len(problems))
        self.assertIn('missing 1 genome(s)', problems[0])

    # the cheap mode

    def testLengthsOnlyCatchesTruncation(self):
        cut = [(n, s[:-1]) if n == 'chr2' else (n, s) for n, s in self.records]
        problems = self.check(cut, lengths_only=True)
        self.assertEqual(1, len(problems))
        self.assertIn('changed the length', problems[0])

    def testLengthsOnlyIgnoresSubstitution(self):
        def substitute(sequence):
            return sequence[:10] + ('T' if sequence[10] != 'T' else 'A') + sequence[11:]
        changed = [(n, substitute(s)) for n, s in self.records]
        self.assertEqual([], self.check(changed, lengths_only=True))


if __name__ == '__main__':
    unittest.main()


class _FakeStore:
    """Just enough of a toil file store for validate_hal_export()."""
    def __init__(self):
        self.log = []

    def readGlobalFile(self, file_id, path, symlink=False):
        shutil.copyfile(file_id, path)

    def deleteLocalFile(self, file_id):
        pass

    def logToMaster(self, message):
        self.log.append(message)


class _FakeJob:
    def __init__(self):
        self.fileStore = _FakeStore()
        self.memory = None

    @property
    def log(self):
        return self.fileStore.log


class _FakeTree:
    """A tree whose node ids are the genome names themselves."""
    def __init__(self, names):
        self.names = names

    def breadthFirstTraversal(self, root=None):
        return list(self.names)

    def getName(self, node):
        return node
