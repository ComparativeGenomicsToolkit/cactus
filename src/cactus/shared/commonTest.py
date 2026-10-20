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
from cactus.shared.common import cactus_call, ChildTreeJob
from cactus.shared.common import makeURL, normalize_pangenome_name, split_pangenome_name, resolve_pangenome_name, pangenome_file_name
from cactus.progressive.seqFile import SeqFile
from cactus.refmap.cactus_minigraph import check_sample_names
from toil.jobStores.fileJobStore import FileJobStore
from urllib.parse import urlparse

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
    def testCactusCall(self):
        inputFile = getTempFile(rootDir=self.tempDir)

        with open("/dev/urandom", "rb") as randText:
            with open(inputFile, 'w') as fh:
                fh.write(b64encode(randText.read(1024)).decode())
        with open(inputFile) as fh:
            input = "".join(fh.read().split("\n"))

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

class TestPangenomeNames(unittest.TestCase):
    """ genome names in a pangenome are PanSN SAMPLE#HAP, with SAMPLE.HAP still read the old way """

    def setUp(self):
        self.tempDir = getTempDirectory(os.getcwd())

    def tearDown(self):
        shutil.rmtree(self.tempDir)

    def testNormalize(self):
        self.assertEqual(normalize_pangenome_name('CHM13'), 'CHM13')
        self.assertEqual(normalize_pangenome_name('HG002#1'), 'HG002#1')
        self.assertEqual(normalize_pangenome_name('HG002#01'), 'HG002#1')
        # the old spelling of a haplotype
        self.assertEqual(normalize_pangenome_name('HG002.2'), 'HG002#2')
        self.assertEqual(normalize_pangenome_name('GCA_054883195.1'), 'GCA_054883195#1')
        self.assertEqual(normalize_pangenome_name('HG002.verkko.1'), 'HG002.verkko#1')
        # any other '.' is part of the sample
        self.assertEqual(normalize_pangenome_name('KOLF2.1J'), 'KOLF2.1J')
        self.assertEqual(normalize_pangenome_name('HG002.verkko'), 'HG002.verkko')
        self.assertEqual(normalize_pangenome_name('GCA_054883195.1#0'), 'GCA_054883195.1#0')
        for bad in ['', '.1', '.hidden', '#1', 'HG002#', 'HG002#x', 'HG002#1#2', 'HG002#1.5', 'a|b', 'a b', 'HG002#¹']:
            with self.assertRaises(RuntimeError, msg=bad):
                normalize_pangenome_name(bad)

    def testSplitAndResolve(self):
        self.assertEqual(split_pangenome_name('CHM13'), ('CHM13', 0))
        self.assertEqual(split_pangenome_name('HG002#2'), ('HG002', 2))
        self.assertEqual(split_pangenome_name('HG002.2'), ('HG002', 2))
        self.assertEqual(split_pangenome_name('GCA_054883195.1#0'), ('GCA_054883195.1', 0))
        names = {'CHM13#0', 'HG002#1', 'GRCh38', 'GCA_054883195.1#0'}
        self.assertEqual(resolve_pangenome_name('HG002#1', names), 'HG002#1')
        self.assertEqual(resolve_pangenome_name('HG002.1', names), 'HG002#1')
        self.assertEqual(resolve_pangenome_name('CHM13', names), 'CHM13#0')
        self.assertEqual(resolve_pangenome_name('GRCh38#0', names), 'GRCh38')
        self.assertEqual(resolve_pangenome_name('GCA_054883195.1#0', names), 'GCA_054883195.1#0')
        # that's haplotype 1 of GCA_054883195, which isn't there
        self.assertIsNone(resolve_pangenome_name('GCA_054883195.1', names))
        self.assertIsNone(resolve_pangenome_name('HG002#2', names))

    def testFileName(self):
        self.assertEqual(pangenome_file_name('HG002#1'), 'HG002.1')
        self.assertEqual(pangenome_file_name('KOLF2.1J'), 'KOLF2.1J')
        self.assertEqual(pangenome_file_name('GCA_054883195.1#0'), 'GCA_054883195.1.0')

    def testMakeURL(self):
        # Toil reads a raw '#' as the start of a fragment, so makeURL quotes it, and Toil unquotes it
        for path in ['/out/HG002#1.fa', '/out/a b%20c.fa', '/out/plain.fa']:
            url = makeURL(path)
            self.assertNotIn('#', url)
            self.assertEqual(FileJobStore._extract_path_from_url(urlparse(url)), path)
        self.assertEqual(makeURL('s3://bucket/x.fa'), 's3://bucket/x.fa')

    def writeSeqFile(self, lines):
        path = os.path.join(self.tempDir, 'seqfile.txt')
        with open(path, 'w') as seq_file:
            seq_file.write('\n'.join(lines) + '\n')
        return path

    def testSeqFile(self):
        path = self.writeSeqFile(['# a comment', 'CHM13\tchm13.fa', 'HG002.1\tpat.fa', 'HG002#2\tmat.fa',
                                  'KOLF2.1J\tkolf.fa', 'GCA_054883195.1#0\tgca.fa'])
        # progressive names are left alone
        self.assertEqual(SeqFile(path).seqOrder, ['CHM13', 'HG002.1', 'HG002#2', 'KOLF2.1J', 'GCA_054883195.1#0'])
        seq_file = SeqFile(path, pangenome=True)
        expected = ['CHM13', 'HG002#1', 'HG002#2', 'KOLF2.1J', 'GCA_054883195.1#0']
        self.assertEqual(seq_file.seqOrder, expected)
        self.assertEqual(seq_file.pathMap['HG002#1'], 'pat.fa')
        self.assertEqual(sorted(seq_file.tree.getName(leaf) for leaf in seq_file.tree.getLeaves()), sorted(expected))
        # and the new names survive being written back out and read again
        path2 = os.path.join(self.tempDir, 'seqfile2.txt')
        with open(path2, 'w') as seq_file2:
            seq_file2.write(str(seq_file))
        self.assertEqual(sorted(SeqFile(path2, pangenome=True).pathMap), sorted(expected))

    def testSeqFileTree(self):
        path = self.writeSeqFile(['(CHM13:1,(HG002.1:1,HG002.2:1)anc:1);', 'CHM13 chm13.fa', 'HG002.1 pat.fa', 'HG002.2 mat.fa'])
        seq_file = SeqFile(path, pangenome=True)
        self.assertEqual(sorted(seq_file.tree.getName(leaf) for leaf in seq_file.tree.getLeaves()), ['CHM13', 'HG002#1', 'HG002#2'])

    def testSeqFileDuplicates(self):
        for names in [('HG002.1', 'HG002#1'), ('CHM13', 'CHM13#0')]:
            path = self.writeSeqFile(['{}\ta.fa'.format(names[0]), '{}\tb.fa'.format(names[1]), 'X\tx.fa'])
            with self.assertRaises(RuntimeError, msg=str(names)):
                SeqFile(path, pangenome=True)

    def testCheckSampleNames(self):
        names = ['CHM13', 'HG002#1', 'HG002#2', 'GCA_054883195.1#0']
        self.assertEqual(check_sample_names(names, ['CHM13']), ['CHM13'])
        self.assertEqual(check_sample_names(names, 'CHM13#0'), 'CHM13')
        self.assertEqual(check_sample_names(names, ['GCA_054883195.1#0', 'CHM13']), ['GCA_054883195.1#0', 'CHM13'])
        with self.assertRaises(RuntimeError):
            check_sample_names(names, ['GRCh38'])
        # only the first reference has to be there (a chromosome's seqfile can leave the others out)
        self.assertEqual(check_sample_names(names, ['CHM13', 'GRCh38.0']), ['CHM13', 'GRCh38#0'])
        # a reference is a whole vg sample, so the other haplotype of its sample can't be there
        with self.assertRaises(RuntimeError):
            check_sample_names(names, ['HG002.1'])
        self.assertEqual(check_sample_names(['CHM13', 'HG002#1'], ['HG002.1']), ['HG002#1'])
        # nor a sample it prefixes
        with self.assertRaises(RuntimeError):
            check_sample_names(['CHM13', 'CHM13v2'], ['CHM13'])
        # dots in the reference's sample don't prefix anything they shouldn't
        self.assertEqual(check_sample_names(['KOLF2.1J', 'KOLF2.2J', 'KOLF22'], 'KOLF2.1J'), 'KOLF2.1J')

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
