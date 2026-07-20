#!/usr/bin/env python3

"""
Unit tests for the pure (non-Toil) logic of cactus-panpatch: the fasta-header sanitizer that has to
match cactus_sanitizeFastaHeaders, the sample-name helpers, and the run resolution that turns a
seqfile (or chromfile) plus options into the list of cactus-pangenome+panpatch runs.

These are fast and offline: make_runs only reads the seqfile text (SeqFile does not open the fastas),
so the fasta paths below are dummies.  The end-to-end pipeline is covered by evolverTest.py.
"""

import os
import unittest
import tempfile
import shutil
from argparse import Namespace

from cactus.refmap.cactus_panpatch import sanitize_contig_name, sample_base, sample_hap
from cactus.refmap.cactus_panpatch import make_runs, panpatch_validate_options

GRAPH_EVENT = '_MINIGRAPH_'

class TestPanpatchUnit(unittest.TestCase):
    def setUp(self):
        self.tempDir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tempDir, ignore_errors=True)

    def _seqfile(self, name, rows):
        """ write a seqfile ("<sample> <fasta>" per row) to tempDir and return its path.  the fasta
        paths are dummies -- make_runs never opens them """
        path = os.path.join(self.tempDir, name)
        with open(path, 'w') as f:
            for sample, fasta in rows:
                f.write('{}\t{}\n'.format(sample, fasta))
        return path

    def _options(self, seqFile, reference=None, sample=None, outName=None, batch=False,
                 noSplit=False, defaultSample=None):
        """ a minimal options namespace with just what panpatch_validate_options / make_runs read.
        panpatch_validate_options fills in referenceFree """
        opts = Namespace(seqFile=seqFile, reference=reference, sample=sample, outName=outName,
                         panpatchBatch=batch, noSplit=noSplit, defaultSample=defaultSample)
        panpatch_validate_options(opts)
        return opts

    def _by_name(self, runs):
        return {run['name']: run for run in runs}

    # ---- sanitize_contig_name (mirror of cactus_sanitizeFastaHeaders -p) ----

    def test_sanitize_contig_name(self):
        # expected values were confirmed against the cactus_sanitizeFastaHeaders -p binary
        cases = {
            'chr1'                    : 'chr1',                 # no-op
            'chr1 some description'   : 'chr1',                 # cut at whitespace
            'GRCh38#0#chr1'           : 'chr1',                 # PanSN -> CONTIG field
            'HG002#2#chr1#extra'      : 'extra',               # only the last '#' field survives
            'chr1:1000-2000'          : 'chr1_sub_999_2000',   # :start-end -> _sub_(start-1)_end
            'HG002#1#chr1:1000-2000'  : 'chr1_sub_999_2000',   # PanSN strip then range
            'chr1:0-500'              : 'chr1_sub_0_500',       # start 0 stays 0 (not -1)
            'scaffold:foo'            : 'scaffold_foo',         # non-range ':' -> '_'
            'a:b:100-200'             : 'a_b_sub_99_200',       # earlier ':' replaced, range converted
            'chr1:2000-1000'          : 'chr1_2000-1000',       # end < start: not a range, ':' -> '_'
            'chr1:10a-20'             : 'chr1_10a-20',          # non-numeric start: not a range
        }
        for header, expected in cases.items():
            self.assertEqual(sanitize_contig_name(header), expected,
                             'sanitize_contig_name({!r})'.format(header))

    # ---- sample_base / sample_hap (only a numeric extension is a haplotype) ----

    def test_sample_base_and_hap(self):
        self.assertEqual(sample_base('HG002.verkko.1'), 'HG002.verkko')  # strip numeric suffix only
        self.assertEqual(sample_base('HG002.verkko'), 'HG002.verkko')    # dotted base is kept
        self.assertEqual(sample_base('S288C'), 'S288C')
        self.assertEqual(sample_hap('HG002.verkko.1'), 1)
        self.assertEqual(sample_hap('HG002.verkko.2'), 2)
        self.assertIsNone(sample_hap('HG002.verkko'))                    # haploid: no suffix
        self.assertEqual(sample_hap('CHM13.0'), 0)

    # ---- make_runs: reference-free (default) ----

    def test_make_runs_reference_free(self):
        # a diploid target whose base name itself contains dots (HG002.verkko), so the .N rename has
        # to flatten every dot: HG002.verkko.1 -> HG002_verkko_1
        seqfile = self._seqfile('rf.seqfile', [
            ('HG002.verkko.1', 'v1.fa'), ('HG002.verkko.2', 'v2.fa'),
            ('HG002.hifiasm.1', 'h1.fa'), ('HG002.hifiasm.2', 'h2.fa')])
        options = self._options(seqfile)   # no --reference -> reference-free
        self.assertTrue(options.referenceFree)

        runs = self._by_name(make_runs(options, GRAPH_EVENT))
        # one run per target haplotype
        self.assertEqual(set(runs), {'HG002.verkko.hap1', 'HG002.verkko.hap2'})

        run1 = runs['HG002.verkko.hap1']
        self.assertTrue(run1['referenceFree'])
        self.assertEqual(run1['hap'], 1)
        self.assertEqual(run1['ploidy'], 1)
        # the target haplotype is renamed to a suffix-free (haploid) sample and is its own reference
        self.assertEqual(run1['panpatchReference'], 'HG002_verkko_1')
        self.assertEqual(run1['reference'], ['HG002_verkko_1'])
        # target first (renamed), then the (un-renamed) donor
        self.assertEqual(run1['samples'], ['HG002_verkko_1', 'HG002.hifiasm'])
        # the graph is the renamed target haplotype plus both donor haplotypes
        self.assertEqual(run1['seq_rows'],
                         [('HG002_verkko_1', 'v1.fa'), ('HG002.hifiasm.1', 'h1.fa'), ('HG002.hifiasm.2', 'h2.fa')])
        # panpatch writes hap0 for the renamed haploid target; target_haps keys the rescue fasta by it
        self.assertEqual(run1['target_haps'], [(0, 'v1.fa')])

        run2 = runs['HG002.verkko.hap2']
        self.assertEqual(run2['panpatchReference'], 'HG002_verkko_2')
        self.assertEqual(run2['seq_rows'][0], ('HG002_verkko_2', 'v2.fa'))
        self.assertEqual(run2['target_haps'], [(0, 'v2.fa')])

    # ---- make_runs: external reference ----

    def test_make_runs_reference(self):
        seqfile = self._seqfile('ref.seqfile', [
            ('CHM13', 'chm13.fa'),
            ('HG002.verkko.1', 'v1.fa'), ('HG002.verkko.2', 'v2.fa'),
            ('HG002.hifiasm.1', 'h1.fa'), ('HG002.hifiasm.2', 'h2.fa')])
        options = self._options(seqfile, reference=['CHM13'])
        self.assertFalse(options.referenceFree)

        runs = make_runs(options, GRAPH_EVENT)
        self.assertEqual(len(runs), 1)   # not per-haplotype
        run = runs[0]
        self.assertFalse(run['referenceFree'])
        self.assertEqual(run['name'], 'HG002.verkko')
        self.assertEqual(run['ploidy'], 2)
        # cactus takes the reference as named in the seqfile; panpatch wants the graph (base) name
        self.assertEqual(run['reference'], ['CHM13'])
        self.assertEqual(run['panpatchReference'], 'CHM13')
        # the reference is not a patch candidate; target first, then donor (base names)
        self.assertEqual(run['samples'], ['HG002.verkko', 'HG002.hifiasm'])
        # reference is placed first in the graph, as cactus-pangenome wants
        self.assertEqual(run['seq_rows'][0], ('CHM13', 'chm13.fa'))
        # both target haplotypes are kept, keyed by their real PanSN haplotype (.1 -> 1, .2 -> 2)
        self.assertEqual(sorted(run['target_haps']), [(1, 'v1.fa'), (2, 'v2.fa')])

    def test_make_runs_reference_and_sample_override(self):
        # --sample chooses the target (and donor order) explicitly, overriding seqfile order
        seqfile = self._seqfile('ref2.seqfile', [
            ('CHM13', 'chm13.fa'),
            ('A.1', 'a1.fa'), ('A.2', 'a2.fa'),
            ('B.1', 'b1.fa'), ('B.2', 'b2.fa')])
        options = self._options(seqfile, reference=['CHM13'], sample=['B', 'A'])
        run = make_runs(options, GRAPH_EVENT)[0]
        self.assertEqual(run['name'], 'B')
        self.assertEqual(run['samples'], ['B', 'A'])   # B patched, A as donor

    # ---- make_runs: --batch ----

    def test_make_runs_batch(self):
        sf1 = self._seqfile('s1.seqfile', [
            ('S1.verkko.1', 'v1.fa'), ('S1.verkko.2', 'v2.fa'),
            ('S1.hifi.1', 'h1.fa'), ('S1.hifi.2', 'h2.fa')])
        sf2 = self._seqfile('s2.seqfile', [
            ('S2.verkko.1', 'v1.fa'), ('S2.verkko.2', 'v2.fa'),
            ('S2.hifi.1', 'h1.fa'), ('S2.hifi.2', 'h2.fa')])
        chromfile = self._seqfile('chromfile.txt', [('SAMP1', sf1), ('SAMP2', sf2)])
        options = self._options(chromfile, batch=True)

        runs = self._by_name(make_runs(options, GRAPH_EVENT))
        # two samples x two haplotypes, named from the chromfile
        self.assertEqual(set(runs), {'SAMP1.hap1', 'SAMP1.hap2', 'SAMP2.hap1', 'SAMP2.hap2'})
        self.assertEqual(runs['SAMP1.hap1']['panpatchReference'], 'S1_verkko_1')
        self.assertEqual(runs['SAMP2.hap2']['panpatchReference'], 'S2_verkko_2')

    # ---- validation rejections ----

    def test_reject_reference_with_haplotype_suffix(self):
        seqfile = self._seqfile('r.seqfile', [('CHM13.1', 'c.fa'), ('A.1', 'a.fa')])
        for bad in (['CHM13.1'], ['CHM13.2'], ['OK', 'CHM13.1']):
            with self.assertRaises(RuntimeError):
                self._options(seqfile, reference=bad)
        # a .0 (or suffix-less) reference is allowed
        self._options(seqfile, reference=['CHM13.0'])

    def test_reject_nosplit(self):
        seqfile = self._seqfile('n.seqfile', [('A.1', 'a.fa'), ('B.1', 'b.fa')])
        with self.assertRaises(RuntimeError):
            self._options(seqfile, noSplit=True)

    def test_reject_batch_with_sample_or_outname(self):
        seqfile = self._seqfile('b.seqfile', [('A.1', 'a.fa'), ('B.1', 'b.fa')])
        with self.assertRaises(RuntimeError):
            self._options(seqfile, batch=True, sample=['A'])
        with self.assertRaises(RuntimeError):
            self._options(seqfile, batch=True, outName='foo')

    def test_reject_missing_seqfile(self):
        with self.assertRaises(RuntimeError):
            self._options(None)

    def test_reject_reference_not_in_seqfile(self):
        seqfile = self._seqfile('m.seqfile', [('A.1', 'a.fa'), ('B.1', 'b.fa')])
        options = self._options(seqfile, reference=['CHM13'])
        with self.assertRaises(RuntimeError):
            make_runs(options, GRAPH_EVENT)

    def test_reject_sample_not_in_seqfile(self):
        seqfile = self._seqfile('s.seqfile', [('A.1', 'a.fa'), ('B.1', 'b.fa')])
        options = self._options(seqfile, sample=['NOPE'])
        with self.assertRaises(RuntimeError):
            make_runs(options, GRAPH_EVENT)

if __name__ == '__main__':
    unittest.main()
