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
from cactus.shared.common import pangenome_seqfile_name
from cactus.progressive.seqFile import SeqFile
from cactus.refmap.cactus_minigraph import check_sample_names
from cactus.refmap.cactus_graphmap_join import graphmap_join_resolve_references
from toil.jobStores.fileJobStore import FileJobStore
from urllib.parse import urlparse
from argparse import Namespace

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
    """ genome names in a pangenome are PanSN SAMPLE#HAP, with SAMPLE.HAP still read the old way in a
    seqfile (dot_haplotypes) and nowhere else """

    def setUp(self):
        self.tempDir = getTempDirectory(os.getcwd())

    def tearDown(self):
        shutil.rmtree(self.tempDir)

    def testNormalize(self):
        for dots in (False, True):
            self.assertEqual(normalize_pangenome_name('CHM13', dots), 'CHM13')
            self.assertEqual(normalize_pangenome_name('HG002#1', dots), 'HG002#1')
            self.assertEqual(normalize_pangenome_name('HG002#01', dots), 'HG002#1')
            # haplotype 0 is spelled without its #0
            self.assertEqual(normalize_pangenome_name('CHM13#0', dots), 'CHM13')
            self.assertEqual(normalize_pangenome_name('GCA_054883195.1#0', dots), 'GCA_054883195.1')
            # any other '.' is part of the sample
            self.assertEqual(normalize_pangenome_name('KOLF2.1J', dots), 'KOLF2.1J')
            self.assertEqual(normalize_pangenome_name('HG002.verkko', dots), 'HG002.verkko')
            for bad in ['', '.1', '.hidden', '#1', 'HG002#', 'HG002#x', 'HG002#1#2', 'HG002#1.5', 'a|b', 'a b', 'HG002#¹']:
                with self.assertRaises(RuntimeError, msg=bad):
                    normalize_pangenome_name(bad, dots)
        # the old spelling of a haplotype, which is only read that way with dot_haplotypes
        for name, pansn in [('HG002.2', 'HG002#2'), ('GCA_054883195.1', 'GCA_054883195#1'),
                            ('HG002.verkko.1', 'HG002.verkko#1'), ('HG002.0', 'HG002')]:
            self.assertEqual(normalize_pangenome_name(name, dot_haplotypes=True), pansn)
            self.assertEqual(normalize_pangenome_name(name), name)

    def testSplitAndResolve(self):
        self.assertEqual(split_pangenome_name('CHM13'), ('CHM13', 0))
        self.assertEqual(split_pangenome_name('HG002#2'), ('HG002', 2))
        # past the seqfile, a '.' is never a haplotype
        self.assertEqual(split_pangenome_name('HG002.2'), ('HG002.2', 0))
        self.assertEqual(split_pangenome_name('GCA_054883195.1'), ('GCA_054883195.1', 0))
        with self.assertRaises(RuntimeError):
            split_pangenome_name('HG002#x')
        names = {'CHM13', 'HG002#1', 'GRCh38', 'GCA_054883195.1'}
        for dots in (False, True):
            self.assertEqual(resolve_pangenome_name('HG002#1', names, dots), 'HG002#1')
            self.assertEqual(resolve_pangenome_name('CHM13#0', names, dots), 'CHM13')
            # a genome that has the name is found by it, whatever dot_haplotypes would make of it
            self.assertEqual(resolve_pangenome_name('GCA_054883195.1', names, dots), 'GCA_054883195.1')
            self.assertEqual(resolve_pangenome_name('GCA_054883195.1#0', names, dots), 'GCA_054883195.1')
            self.assertIsNone(resolve_pangenome_name('HG002#2', names, dots))
        # otherwise the name is read as the seqfile's names were
        self.assertEqual(resolve_pangenome_name('HG002.1', names, True), 'HG002#1')
        self.assertIsNone(resolve_pangenome_name('HG002.1', names, False))

    def testFileName(self):
        self.assertEqual(pangenome_file_name('HG002#1'), 'HG002.1')
        self.assertEqual(pangenome_file_name('KOLF2.1J'), 'KOLF2.1J')
        self.assertEqual(pangenome_file_name('GCA_054883195.1'), 'GCA_054883195.1')

    def testSeqFileName(self):
        # a name cactus writes into a seqfile reads back as the same genome whatever dot_haplotypes
        for name in ['CHM13', 'HG002#1', 'KOLF2.1J', 'HG002.verkko', 'GCA_054883195.1']:
            for dots in (False, True):
                self.assertEqual(normalize_pangenome_name(pangenome_seqfile_name(name), dots), name)
        self.assertEqual(pangenome_seqfile_name('GCA_054883195.1'), 'GCA_054883195.1#0')
        self.assertEqual(pangenome_seqfile_name('KOLF2.1J'), 'KOLF2.1J')

    def testMakeURL(self):
        # Toil reads a raw '#' as the start of a fragment, so makeURL quotes it, and Toil unquotes it
        for path in ['/out/HG002#1.fa', '/out/a b%20c.fa', '/out/plain.fa']:
            url = makeURL(path)
            self.assertNotIn('#', url)
            self.assertEqual(FileJobStore._extract_path_from_url(urlparse(url)), path)
        self.assertEqual(makeURL('s3://bucket/x.fa'), 's3://bucket/x.fa')

    def writeSeqFile(self, lines, name='seqfile.txt'):
        path = os.path.join(self.tempDir, name)
        with open(path, 'w') as seq_file:
            seq_file.write('\n'.join(lines) + '\n')
        return path

    def testSeqFile(self):
        path = self.writeSeqFile(['# a comment', 'CHM13\tchm13.fa', 'HG002.1\tpat.fa', 'HG002#2\tmat.fa',
                                  'KOLF2.1J\tkolf.fa', 'GCA_054883195.1#0\tgca.fa'])
        # progressive names are left alone
        self.assertEqual(SeqFile(path).seqOrder, ['CHM13', 'HG002.1', 'HG002#2', 'KOLF2.1J', 'GCA_054883195.1#0'])
        seq_file = SeqFile(path, pangenome=True, dot_haplotypes=True)
        expected = ['CHM13', 'HG002#1', 'HG002#2', 'KOLF2.1J', 'GCA_054883195.1']
        self.assertEqual(seq_file.seqOrder, expected)
        self.assertEqual(seq_file.pathMap['HG002#1'], 'pat.fa')
        self.assertEqual(sorted(seq_file.tree.getName(leaf) for leaf in seq_file.tree.getLeaves()), sorted(expected))
        # without dot_haplotypes, HG002.1 is a sample of its own
        self.assertEqual(SeqFile(path, pangenome=True, dot_haplotypes=False).seqOrder,
                         ['CHM13', 'HG002.1', 'HG002#2', 'KOLF2.1J', 'GCA_054883195.1'])
        # and the new names survive being written back out the way cactus writes them (cactus-graphmap-split's
        # seqfiles), and read again whatever dot_haplotypes
        path2 = self.writeSeqFile(['{}\t{}'.format(pangenome_seqfile_name(name), seq_file.pathMap[name]) for name in expected],
                                  'seqfile2.txt')
        for dots in (False, True):
            self.assertEqual(SeqFile(path2, pangenome=True, dot_haplotypes=dots).seqOrder, expected)
        # a pangenome seqfile can't be read without saying which way
        with self.assertRaises(AssertionError):
            SeqFile(path, pangenome=True)
        # the old spelling is called out either way: rewritten with dot_haplotypes, and without it, where it
        # looks like haplotypes of one sample, taken for samples of their own
        with self.assertLogs('toil.statsAndLogging', level='WARNING') as logs:
            SeqFile(path, pangenome=True, dot_haplotypes=True)
        self.assertIn('HG002.1 -> HG002#1', logs.output[0])
        path3 = self.writeSeqFile(['CHM13\tchm13.fa', 'HG002.1\tpat.fa', 'HG002.2\tmat.fa', 'GCA_054883195.1\tgca.fa'],
                                  'seqfile3.txt')
        with self.assertLogs('toil.statsAndLogging', level='WARNING') as logs:
            SeqFile(path3, pangenome=True, dot_haplotypes=False)
        self.assertEqual(len(logs.output), 1)
        self.assertIn('HG002.1 HG002.2 are each a haploid sample', logs.output[0])

    def testSeqFileTree(self):
        path = self.writeSeqFile(['(CHM13:1,(HG002.1:1,HG002.2:1)anc:1);', 'CHM13 chm13.fa', 'HG002.1 pat.fa', 'HG002.2 mat.fa'])
        seq_file = SeqFile(path, pangenome=True, dot_haplotypes=True)
        self.assertEqual(sorted(seq_file.tree.getName(leaf) for leaf in seq_file.tree.getLeaves()), ['CHM13', 'HG002#1', 'HG002#2'])
        seq_file = SeqFile(path, pangenome=True, dot_haplotypes=False)
        self.assertEqual(sorted(seq_file.tree.getName(leaf) for leaf in seq_file.tree.getLeaves()), ['CHM13', 'HG002.1', 'HG002.2'])

    def testSeqFileDuplicates(self):
        # the same genome twice, and two genomes that would give the same file names
        for names, dots in [(('HG002.1', 'HG002#1'), True), (('CHM13', 'CHM13#0'), True), (('CHM13', 'CHM13#0'), False),
                            (('GCA_054883195.1#0', 'GCA_054883195#1'), True), (('HG002.1', 'HG002#1'), False)]:
            path = self.writeSeqFile(['{}\ta.fa'.format(names[0]), '{}\tb.fa'.format(names[1]), 'X\tx.fa'])
            with self.assertRaises(RuntimeError, msg=str(names)):
                SeqFile(path, pangenome=True, dot_haplotypes=dots)

    def testCheckSampleNames(self):
        names = ['CHM13', 'HG002#1', 'HG002#2', 'GCA_054883195.1']
        for dots in (False, True):
            self.assertEqual(check_sample_names(names, ['CHM13'], dots), ['CHM13'])
            self.assertEqual(check_sample_names(names, 'CHM13#0', dots), 'CHM13')
            # a haplotype 0 genome whose name ends in .N is found by its name, or with its #0
            self.assertEqual(check_sample_names(names, ['GCA_054883195.1', 'CHM13'], dots), ['GCA_054883195.1', 'CHM13'])
            self.assertEqual(check_sample_names(names, ['GCA_054883195.1#0', 'CHM13'], dots), ['GCA_054883195.1', 'CHM13'])
            with self.assertRaises(RuntimeError):
                check_sample_names(names, ['GRCh38'], dots)
            # only the first reference has to be there (a chromosome's seqfile can leave the others out), and
            # one that isn't keeps its name (bar a #0)
            self.assertEqual(check_sample_names(names, ['CHM13', 'GRCh38#0'], dots), ['CHM13', 'GRCh38'])
            self.assertEqual(check_sample_names(names, ['CHM13', 'GRCh38.1'], dots), ['CHM13', 'GRCh38.1'])
            # a reference is a whole vg sample, so the other haplotype of its sample can't be there
            with self.assertRaises(RuntimeError):
                check_sample_names(names, ['HG002#1'], dots)
            # nor a sample it prefixes
            with self.assertRaises(RuntimeError):
                check_sample_names(['CHM13', 'CHM13v2'], ['CHM13'], dots)
            # dots in the reference's sample don't prefix anything they shouldn't
            self.assertEqual(check_sample_names(['KOLF2.1J', 'KOLF2.2J', 'KOLF22'], 'KOLF2.1J', dots), 'KOLF2.1J')
        # the old spelling of a haplotype finds it only with dot_haplotypes, as in the seqfile
        self.assertEqual(check_sample_names(['CHM13', 'HG002#1'], ['HG002.1'], True), ['HG002#1'])
        with self.assertRaises(RuntimeError):
            check_sample_names(['CHM13', 'HG002#1'], ['HG002.1'], False)
        # once resolved against one seqfile, a reference comes through another one unchanged, even one that
        # leaves it out (cactus-pangenome, and --batch, resolve against the whole seqfile first)
        resolved = check_sample_names(names, ['CHM13', 'GCA_054883195.1'], True)
        self.assertEqual(check_sample_names(['CHM13', 'X'], resolved, True), ['CHM13', 'GCA_054883195.1'])
        self.assertEqual(check_sample_names(names, resolved, True), ['CHM13', 'GCA_054883195.1'])

    def testJoinReferences(self):
        def resolve(reference, vcf_reference=None, dots=True, resolved=False, bypass=False):
            options = Namespace(reference=list(reference), vcfReference=list(vcf_reference) if vcf_reference else None,
                                bypass=bypass, vcf=['clip'])
            graphmap_join_resolve_references(options, dots, references_resolved=resolved)
            return options.reference, options.vcfReference
        # cactus-graphmap-join has no seqfile, so it reads the references the way a seqfile's names are read
        self.assertEqual(resolve(['CHM13#0', 'HG002.1']), (['CHM13', 'HG002#1'], ['CHM13']))
        self.assertEqual(resolve(['CHM13#0', 'HG002.1'], dots=False), (['CHM13', 'HG002.1'], ['CHM13']))
        self.assertEqual(resolve(['GCA_054883195.1#0']), (['GCA_054883195.1'], ['GCA_054883195.1']))
        # and says when that makes a reference a haplotype, as it can't tell the two apart
        with self.assertLogs('toil.statsAndLogging', level='WARNING') as logs:
            resolve(['CHM13', 'GCA_054883195.1'])
        self.assertIn('give it as GCA_054883195.1#0', logs.output[0])
        # a --vcfReference is one of the references, found by its name first
        self.assertEqual(resolve(['CHM13', 'GCA_054883195.1#0'], ['GCA_054883195.1']),
                         (['CHM13', 'GCA_054883195.1'], ['GCA_054883195.1']))
        self.assertEqual(resolve(['CHM13', 'HG002.1'], ['HG002.1']), (['CHM13', 'HG002#1'], ['HG002#1']))
        with self.assertRaises(RuntimeError):
            resolve(['CHM13'], ['GRCh38'])
        # references cactus-pangenome has resolved against its seqfile are left alone
        self.assertEqual(resolve(['GCA_054883195.1', 'CHM13'], resolved=True), (['GCA_054883195.1', 'CHM13'], ['GCA_054883195.1']))
        # in bypass mode, only the first reference can have a VCF
        with self.assertRaises(RuntimeError):
            resolve(['CHM13', 'GRCh38'], ['GRCh38'], bypass=True)

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
