#!/usr/bin/env python3

"""
Unit tests for the pure (non-Toil) logic of cactus-panpatch: the fasta-header sanitizer that has to
match cactus_sanitizeFastaHeaders, the sample-name helpers, and the run resolution that turns a
seqfile (or chromfile) plus options into the list of cactus-pangenome+panpatch runs.

These are fast and offline: make_runs only reads the seqfile text (SeqFile does not open the fastas),
so the fasta paths below are dummies.  The end-to-end pipeline is covered by evolverTest.py.
"""

import os
import re
import gzip
import random
import unittest
import tempfile
import shutil
from argparse import Namespace

from cactus.refmap.cactus_panpatch import sanitize_contig_name, sample_base, sample_hap
from cactus.refmap.cactus_panpatch import make_runs, panpatch_validate_options
from cactus.refmap.cactus_panpatch import read_error_bed_manifest, parse_error_bed, intervals_overlap
from cactus.refmap.cactus_panpatch import parse_bed_blocks, make_panpatch_error_bed
from cactus.refmap.cactus_panpatch import tag_report_gap_origin, summarize_report, write_batch_summary, concat_reports
from cactus.refmap.cactus_panpatch import contig_is_t2t
from cactus.refmap.cactus_panpatch import target_error_intervals_used, panpatch_option_value
from cactus.refmap.cactus_panpatch import clean_path_name, clean_report_names, clean_bed_names, run_renames
from cactus.refmap.cactus_panpatch import output_record_name, rename_fasta_records, missing_target_contigs
from cactus.refmap.cactus_panpatch import rescue_dropped_contigs, PANPATCH_REPORT_COLUMNS
from cactus.refmap.cactus_panpatch import target_error_bed_keys, error_interval_sequences, covered_fraction
from cactus.refmap.cactus_panpatch import tandem_filter_bed, filter_target_error_beds

GRAPH_EVENT = '_MINIGRAPH_'
# <graphmap dotHaplotypes> as the default config has it, which is what the seqfiles below, with their
# old SAMPLE.N names, are written for
DOT_HAPLOTYPES = True

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
            'id=HG002#1|chr1'         : 'chr1',                 # already sanitized: '#' in the genome
        }
        for header, expected in cases.items():
            self.assertEqual(sanitize_contig_name(header), expected,
                             'sanitize_contig_name({!r})'.format(header))

    # ---- sample_base / sample_hap (only a #N is a haplotype) ----

    def test_sample_base_and_hap(self):
        # PanSN: '#' gives the haplotype, and a '.' is part of the sample name.  these take genome names as
        # the seqfile's are normalized, so the old SAMPLE.N spelling has been rewritten to SAMPLE#N already
        self.assertEqual(sample_base('HG002.verkko#1'), 'HG002.verkko')
        self.assertEqual(sample_hap('HG002.verkko#1'), 1)
        self.assertEqual(sample_hap('HG002.verkko#2'), 2)
        self.assertEqual(sample_base('HG002.verkko'), 'HG002.verkko')    # dotted base is kept
        self.assertIsNone(sample_hap('HG002.verkko'))                    # haploid: no haplotype
        self.assertEqual(sample_base('S288C'), 'S288C')
        self.assertEqual(sample_base('KOLF2.1J'), 'KOLF2.1J')
        self.assertIsNone(sample_hap('KOLF2.1J'))
        self.assertEqual(sample_base('GCA_054883195.1'), 'GCA_054883195.1')
        self.assertIsNone(sample_hap('GCA_054883195.1'))
        self.assertEqual(sample_base('HG002.verkko.1'), 'HG002.verkko.1')
        self.assertIsNone(sample_hap('HG002.verkko.1'))

    # ---- make_runs: reference-free (default) ----

    def test_make_runs_reference_free(self):
        # a diploid target whose base name itself contains dots (HG002.verkko), so the .N rename has
        # to flatten every dot: HG002.verkko.1 -> HG002_verkko_1
        seqfile = self._seqfile('rf.seqfile', [
            ('HG002.verkko.1', 'v1.fa'), ('HG002.verkko.2', 'v2.fa'),
            ('HG002.hifiasm.1', 'h1.fa'), ('HG002.hifiasm.2', 'h2.fa')])
        options = self._options(seqfile)   # no --reference -> reference-free
        self.assertTrue(options.referenceFree)

        runs = self._by_name(make_runs(options, GRAPH_EVENT, DOT_HAPLOTYPES))
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
        # the graph is the renamed target haplotype plus both donor haplotypes (whose old SAMPLE.N
        # names are rewritten to SAMPLE#N as the seqfile is read)
        self.assertEqual(run1['seq_rows'],
                         [('HG002_verkko_1', 'v1.fa'), ('HG002.hifiasm#1', 'h1.fa'), ('HG002.hifiasm#2', 'h2.fa')])
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

        runs = make_runs(options, GRAPH_EVENT, DOT_HAPLOTYPES)
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
        run = make_runs(options, GRAPH_EVENT, DOT_HAPLOTYPES)[0]
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

        runs = self._by_name(make_runs(options, GRAPH_EVENT, DOT_HAPLOTYPES))
        # two samples x two haplotypes, named from the chromfile
        self.assertEqual(set(runs), {'SAMP1.hap1', 'SAMP1.hap2', 'SAMP2.hap1', 'SAMP2.hap2'})
        self.assertEqual(runs['SAMP1.hap1']['panpatchReference'], 'S1_verkko_1')
        self.assertEqual(runs['SAMP2.hap2']['panpatchReference'], 'S2_verkko_2')

    # ---- validation rejections ----

    def test_reject_reference_with_haplotype_suffix(self):
        # CHM13.1 is the old spelling of CHM13#1, one haplotype of a diploid, which can't be the reference
        seqfile = self._seqfile('r.seqfile', [('CHM13.1', 'c.fa'), ('A.1', 'a.fa')])
        for bad in (['CHM13.1'], ['CHM13#1'], ['CHM13.2'], ['OK', 'CHM13.1']):
            with self.assertRaises(RuntimeError):
                make_runs(self._options(seqfile, reference=bad), GRAPH_EVENT, DOT_HAPLOTYPES)
        # with dotHaplotypes="0", it's a haploid sample named CHM13.1
        run = make_runs(self._options(seqfile, reference=['CHM13.1']), GRAPH_EVENT, False)[0]
        self.assertEqual(run['panpatchReference'], 'CHM13.1')
        # a .0, #0 or suffix-less reference is allowed
        seqfile = self._seqfile('r0.seqfile', [('CHM13', 'c.fa'), ('A.1', 'a.fa')])
        for good in (['CHM13'], ['CHM13#0'], ['CHM13.0']):
            run = make_runs(self._options(seqfile, reference=good), GRAPH_EVENT, DOT_HAPLOTYPES)[0]
            self.assertEqual(run['reference'], ['CHM13'])

    def test_make_runs_versioned_reference(self):
        # an accession's version: GCA_054883195.1#0 keeps it with dotHaplotypes, after which the genome is
        # GCA_054883195.1 and found by that name.  with dotHaplotypes="0" it needs no #0 in the first place
        for dots, ref_name in [(True, 'GCA_054883195.1#0'), (False, 'GCA_054883195.1')]:
            seqfile = self._seqfile('v.seqfile', [(ref_name, 'ref.fa'), ('TGT#1', 't1.fa'), ('TGT#2', 't2.fa')])
            for reference in (['GCA_054883195.1'], ['GCA_054883195.1#0']):
                run = make_runs(self._options(seqfile, reference=reference), GRAPH_EVENT, dots)[0]
                self.assertEqual(run['reference'], ['GCA_054883195.1'])
                self.assertEqual(run['panpatchReference'], 'GCA_054883195.1')
                self.assertEqual(run['seq_rows'][0], ('GCA_054883195.1', 'ref.fa'))
                self.assertEqual(sorted(run['target_haps']), [(1, 't1.fa'), (2, 't2.fa')])
        # whereas a seqfile's GCA_054883195.1 is haplotype 1 of GCA_054883195 with dotHaplotypes, which can't be a
        # reference
        seqfile = self._seqfile('v1.seqfile', [('GCA_054883195.1', 'ref.fa'), ('TGT#1', 't1.fa')])
        with self.assertRaises(RuntimeError):
            make_runs(self._options(seqfile, reference=['GCA_054883195.1']), GRAPH_EVENT, True)

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
            make_runs(options, GRAPH_EVENT, DOT_HAPLOTYPES)

    def test_reject_sample_not_in_seqfile(self):
        seqfile = self._seqfile('s.seqfile', [('A.1', 'a.fa'), ('B.1', 'b.fa')])
        options = self._options(seqfile, sample=['NOPE'])
        with self.assertRaises(RuntimeError):
            make_runs(options, GRAPH_EVENT, DOT_HAPLOTYPES)

class _FakeFileStore:
    """ a minimal Toil fileStore backed by local files, so the file-id-taking helpers (mask / revert /
    tag) can be exercised offline.  a "file id" here is just a local path """
    def __init__(self, root):
        self.root = root
        self._n = 0

    def _fresh(self, tag):
        self._n += 1
        return os.path.join(self.root, '{}_{}'.format(tag, self._n))

    def getLocalTempDir(self):
        d = self._fresh('work')
        os.makedirs(d, exist_ok=True)
        return d

    def readGlobalFile(self, file_id, path):
        shutil.copyfile(file_id, path)

    def writeGlobalFile(self, path):
        dst = self._fresh('stored')
        shutil.copyfile(path, dst)
        return dst

    def exportFile(self, file_id, url):
        # makeURL gives file://<abspath> for a local path
        shutil.copyfile(file_id, url[len('file://'):] if url.startswith('file://') else url)

class _FakeJob:
    def __init__(self, root):
        self.fileStore = _FakeFileStore(root)

# panpatch's per-patch report header (mirror of PATCH_TABLE_HEADER in panpatch_main.cpp)
_REPORT_HEADER = ('chrom\thap\ttype\ttarget\ttarget_bp\tdonor\tdonor_bp\treplaced_bp\t'
                  'kmer%\tflankL%\tflankR%\tdecision\treason\ttarget_start\ttarget_end')

def _read_fasta_bytes(path):
    """ {record_name: sequence_bytes} for a (possibly wrapped) fasta """
    recs, name, parts = {}, None, []
    with open(path, 'rb') as f:
        for line in f:
            if line.startswith(b'>'):
                if name is not None:
                    recs[name] = b''.join(parts)
                name, parts = line[1:].rstrip(b'\n').decode(), []
            else:
                parts.append(line.strip())
    if name is not None:
        recs[name] = b''.join(parts)
    return recs

class TestPanpatchErrorBeds(unittest.TestCase):
    """ the assembly-error-BED wrapper logic: manifest parsing, the small pure helpers, the N-mask
    preprocess, the whole-interval target revert, and the gap_origin report tagging """

    def setUp(self):
        self.tempDir = tempfile.mkdtemp()
        self.job = _FakeJob(self.tempDir)

    def tearDown(self):
        shutil.rmtree(self.tempDir, ignore_errors=True)

    def _write(self, name, text, binary=False):
        path = os.path.join(self.tempDir, name)
        with open(path, 'wb' if binary else 'w') as f:
            f.write(text)
        return path

    def _seqfile(self, name, rows):
        path = os.path.join(self.tempDir, name)
        with open(path, 'w') as f:
            for sample, fasta in rows:
                f.write('{}\t{}\n'.format(sample, fasta))
        return path

    def _by_name(self, runs):
        return {run['name']: run for run in runs}

    # ---- pure helpers ----

    def test_intervals_overlap(self):
        ivs = [(10, 20), (30, 40)]
        self.assertTrue(intervals_overlap(15, 16, ivs))
        self.assertTrue(intervals_overlap(5, 11, ivs))         # touches the left edge
        self.assertFalse(intervals_overlap(20, 30, ivs))       # half-open: [20,30) abuts, no overlap
        self.assertFalse(intervals_overlap(0, 10, ivs))
        self.assertFalse(intervals_overlap(0, 5, []))

    def test_parse_error_bed(self):
        bed = self._write('e.bed', 'track name=x\nchr1\t10\t20\nchr1\t30\t40\nGRCh38#0#chr2\t5\t9\n'
                                   'bad\t50\t50\nchr3\t-1\t5\n')  # zero-length and negative dropped
        raw = parse_error_bed(bed)
        self.assertEqual(raw, {'chr1': [(10, 20), (30, 40)], 'GRCh38#0#chr2': [(5, 9)]})
        san = parse_error_bed(bed, sanitize=True)                # PanSN -> CONTIG field
        self.assertEqual(san['chr2'], [(5, 9)])
        self.assertIn('chr1', san)

    def test_read_error_bed_manifest(self):
        a = self._write('a.bed', 'chr1\t0\t5\n')
        b = self._write('b.bed', 'chr2\t0\t5\n')
        man = self._write('manifest.txt', '# a comment\nTGT.1\t{}\n\nDON\t{}\n'.format(a, b))
        got = read_error_bed_manifest(man)
        self.assertEqual(got, {'TGT.1': os.path.abspath(a), 'DON': os.path.abspath(b)})
        # a duplicate event is rejected
        dup = self._write('dup.txt', 'TGT.1\t{}\nTGT.1\t{}\n'.format(a, b))
        with self.assertRaises(RuntimeError):
            read_error_bed_manifest(dup)
        # a missing bed path is rejected
        miss = self._write('miss.txt', 'TGT.1\t{}\n'.format(os.path.join(self.tempDir, 'nope.bed')))
        with self.assertRaises(RuntimeError):
            read_error_bed_manifest(miss)

    def test_parse_bed_blocks(self):
        bed = ('#Patched assembly on chr1 for TGT#1:\n'
               'TGT#1#chr1\t0\t60\t+\n'
               'DON#1#dcontig\t60\t80\t+\n'
               'TGT#1#chr1\t80\t100\t+\n'
               '\n'
               '#Passthrough (unplaced) for TGT#1:\n'
               'TGT#1#unplaced\t0\t50\t+\n'
               '\n'
               '#Patched assembly on chr2 for TGT#2:\n'
               'TGT#2#chr2\t0\t24\t-\n')
        blocks = parse_bed_blocks(self._write('p.bed', bed))
        self.assertEqual(blocks[('chr1', '1')],
                         [('TGT#1#chr1', 0, 60, '+'), ('DON#1#dcontig', 60, 80, '+'),
                          ('TGT#1#chr1', 80, 100, '+')])
        self.assertEqual(blocks[('chr2', '2')], [('TGT#2#chr2', 0, 24, '-')])
        # the passthrough block is not a "#Patched assembly on" block, so it is never captured
        self.assertNotIn(('unplaced', '1'), blocks)
        self.assertEqual(len(blocks), 2)

    # ---- make_panpatch_error_bed ----

    def test_make_panpatch_error_bed(self):
        # the reference-free target is a haploid sample in the graph (TGT_1 -> haplotype 0); a donor keeps
        # its seqfile sample and haplotype, and its contig names (here PanSN-prefixed, as the hprc-r2
        # fastas have them) are sanitized down to the CONTIG field the graph uses
        tgt = self._write('t.bed', 'chr1_haplotype1-0000003\t10\t20\nchr1_haplotype1-0000003\t30\t40\n')
        don = self._write('d.bed', 'HG1#2#CM1.1\t5\t9\n# a comment\nctg:11-20\t0\t4\n')
        out = make_panpatch_error_bed(self.job, {'TGT_1': tgt, 'DON#2': don})
        with open(out) as f:
            lines = f.read().splitlines()
        # assemblies in sample-name order, each BED's intervals in its own order
        self.assertEqual(lines, ['DON#2#CM1.1\t5\t9',
                                 'DON#2#ctg_sub_10_20\t0\t4',
                                 'TGT_1#0#chr1_haplotype1-0000003\t10\t20',
                                 'TGT_1#0#chr1_haplotype1-0000003\t30\t40'])

    def test_make_panpatch_error_bed_empty(self):
        out = make_panpatch_error_bed(self.job, {'TGT_1': self._write('e.bed', '')})
        with open(out) as f:
            self.assertEqual(f.read(), '')

    # ---- tag_report_gap_origin ----

    def test_tag_report_gap_origin(self):
        # a target error bed on hap 1, chr1 (100 kb), covering [30000,50000)
        err = self._write('t1.bed', 'chr1\t30000\t50000\n')
        target_error_bed_ids = {1: err}

        def row(chrom, hap, typ, target, ts, te):
            return '\t'.join([chrom, str(hap), typ, target, '100000', 'DON#1#d', '20', '0',
                              '.', '100.0', '100.0', 'accepted', '.', str(ts), str(te)])
        report = self._write('report.tsv', _REPORT_HEADER + '\n'
                             + row('chr1', 1, 'gap-fill', 'TGT#1#chr1', 35000, 45000) + '\n'   # overlaps -> bed-gap
                             + row('chr1', 1, 'gap-fill', 'TGT#1#chr1', 70000, 80000) + '\n'   # no overlap -> N-gap
                             + row('chr1', 1, 'telomere', 'TGT#1#chr1', 45000, 55000) + '\n'   # overlaps -> bed-gap
                             + row('chr1', 1, 'telomere', 'TGT#1#chr1', 0, 10000) + '\n'       # no overlap (unflagged end) -> .
                             + row('chr1', 2, 'gap-fill', 'TGT#2#chr1', 35000, 45000) + '\n'   # hap 2 has no bed -> N-gap
                             + '#Contig TGT#1#chr1 len=100000\n')

        tag_report_gap_origin(self.job, report, target_error_bed_ids)
        with open(report) as f:
            lines = [l.rstrip('\n') for l in f]

        self.assertEqual(lines[0], _REPORT_HEADER + '\tgap_origin')
        self.assertTrue(lines[1].endswith('\tbed-gap'))
        self.assertTrue(lines[2].endswith('\tN-gap'))
        self.assertTrue(lines[3].endswith('\tbed-gap'))    # telomere over a flagged interval
        self.assertTrue(lines[4].endswith('\t.'))          # telomere over an unflagged end
        self.assertTrue(lines[5].endswith('\tN-gap'))      # hap 2, no error bed
        self.assertEqual(lines[6], '#Contig TGT#1#chr1 len=100000')   # comment passes through, no column

    def test_tag_report_gap_origin_contig_ends(self):
        # contig-end flags, which panpatch ignores on the target, explain no patch: [0,20000) and [20000,25000) join
        # up into one from chr1's start, [95000,100000) runs to its end, and [60000,65000) is an ordinary flag
        err = self._write('t2.bed', 'chr1\t0\t20000\nchr1\t20000\t25000\nchr1\t60000\t65000\nchr1\t95000\t100000\n')

        def row(typ, ts, te):
            return '\t'.join(['chr1', '1', typ, 'TGT#1#chr1', '100000', 'DON#1#d', '20', '0',
                              '.', '100.0', '100.0', 'accepted', '.', str(ts), str(te)])
        rows = [row('telomere', 0, 10000), row('gap-fill', 22000, 23000), row('gap-fill', 61000, 62000),
                row('telomere', 96000, 100000)]
        for slop, want in [(1000, ['.', 'N-gap', 'bed-gap', '.']), (0, ['bed-gap', 'bed-gap', 'bed-gap', 'bed-gap'])]:
            report = self._write('ends{}.tsv'.format(slop), _REPORT_HEADER + '\n' + '\n'.join(rows) + '\n')
            tag_report_gap_origin(self.job, report, {1: err}, slop)
            with open(report) as f:
                got = [l.rstrip('\n').split('\t')[-1] for l in f][1:]
            self.assertEqual(got, want, 'end_slop {}'.format(slop))

    def test_tag_report_gap_origin_long(self):
        # a target interval over max_len (1 Mb by default, as panpatch's --error-bed-max-len) explains no patch
        err = self._write('t3.bed', 'chr1\t100000\t1100001\nchr1\t2000000\t3000000\n')

        def row(ts, te):
            return '\t'.join(['chr1', '1', 'gap-fill', 'TGT#1#chr1', '5000000', 'DON#1#d', '20', '0',
                              '.', '100.0', '100.0', 'accepted', '.', str(ts), str(te)])
        rows = [row(500000, 510000), row(2500000, 2510000)]
        for max_len, want in [(None, ['N-gap', 'bed-gap']), (0, ['bed-gap', 'bed-gap'])]:
            report = self._write('long{}.tsv'.format(max_len), _REPORT_HEADER + '\n' + '\n'.join(rows) + '\n')
            if max_len is None:
                tag_report_gap_origin(self.job, report, {1: err})
            else:
                tag_report_gap_origin(self.job, report, {1: err}, 1000, max_len)
            with open(report) as f:
                got = [l.rstrip('\n').split('\t')[-1] for l in f][1:]
            self.assertEqual(got, want, 'max_len {}'.format(max_len))

    def test_target_error_intervals_used(self):
        # merged first, so a flag split over lines is judged as one; the slop boundary is inclusive, as in panpatch
        self.assertEqual(target_error_intervals_used([(15000, 28999), (0, 15000), (40000, 50000)], 100000, 1000),
                         [(40000, 50000)])
        self.assertEqual(target_error_intervals_used([(1000, 5000), (1001, 5000)], 100000, 1000), [])
        self.assertEqual(target_error_intervals_used([(1001, 5000), (95000, 98999)], 100000, 1000),
                         [(1001, 5000), (95000, 98999)])
        self.assertEqual(target_error_intervals_used([(95000, 99000)], 100000, 1000), [])
        self.assertEqual(target_error_intervals_used([(0, 5000)], 100000, 0), [(0, 5000)])
        # and one longer than max_len (after the merge; the limit itself is kept), unless max_len is 0
        self.assertEqual(target_error_intervals_used([(10000, 20001), (30000, 40000)], 100000, 1000, 10000),
                         [(30000, 40000)])
        self.assertEqual(target_error_intervals_used([(10000, 15000), (15000, 20001)], 100000, 1000, 10000), [])
        self.assertEqual(target_error_intervals_used([(10000, 20001)], 100000, 1000, 0), [(10000, 20001)])

    def test_panpatch_option_value(self):
        self.assertEqual(panpatch_option_value(None, '--error-bed-end-slop', 1000), 1000)
        self.assertEqual(panpatch_option_value('--min-cover 0.9', '--error-bed-end-slop', 1000), 1000)
        self.assertEqual(panpatch_option_value('--error-bed-end-slop 0 --min-cover 0.9', '--error-bed-end-slop', 1000), '0')
        self.assertEqual(panpatch_option_value('--error-bed-end-slop=500', '--error-bed-end-slop', 1000), '500')
        self.assertEqual(panpatch_option_value('--error-bed-end-slop 5 --error-bed-end-slop 7', '--error-bed-end-slop', 1000), '7')
        # as the post-revert telomere measurement reads the threshold: unset -> None, either form -> its value
        self.assertIsNone(panpatch_option_value('--min-cover 0.9', '--telomere-threshold', None))
        self.assertEqual(panpatch_option_value('--min-cover 0.9 --telomere-threshold=0.7', '--telomere-threshold', None), '0.7')

    def test_tag_report_no_beds_still_adds_column(self):
        # with no target error beds, every gap-fill is a genuine N-gap but the column is still added
        def row(typ, ts, te):
            return '\t'.join(['chr1', '1', typ, 'TGT#1#chr1', '100', 'DON#1#d', '20', '0',
                              '.', '100.0', '100.0', 'accepted', '.', str(ts), str(te)])
        report = self._write('r2.tsv', _REPORT_HEADER + '\n' + row('gap-fill', 10, 20) + '\n')
        tag_report_gap_origin(self.job, report, {})
        with open(report) as f:
            lines = [l.rstrip('\n') for l in f]
        self.assertEqual(lines[0], _REPORT_HEADER + '\tgap_origin')
        self.assertTrue(lines[1].endswith('\tN-gap'))

    # ---- --fromGraphs ----

    def test_from_graphs_out_dir(self):
        # --fromGraphs must not write its outputs over the run whose graphs it reads
        seqfile = self._seqfile('g.seqfile', [('TGT.1', 't1.fa'), ('DON.1', 'd1.fa')])
        kept = os.path.join(self.tempDir, 'kept')
        os.makedirs(kept)
        opts = Namespace(seqFile=seqfile, reference=None, sample=None, outName=None, panpatchBatch=False,
                         noSplit=False, defaultSample=None, fromGraphs=kept, outDir=kept + '/')
        with self.assertRaises(RuntimeError):
            panpatch_validate_options(opts)
        opts.outDir = os.path.join(self.tempDir, 'new')
        panpatch_validate_options(opts)

    # ---- make_runs wiring of the error beds ----

    def _options_with_beds(self, seqFile, error_map, reference=None, sample=None):
        opts = Namespace(seqFile=seqFile, reference=reference, sample=sample, outName=None,
                         panpatchBatch=False, noSplit=False, defaultSample=None)
        panpatch_validate_options(opts)
        opts.errorBedMap = error_map
        return opts

    def test_make_runs_error_beds_reference(self):
        seqfile = self._seqfile('r.seqfile', [('CHM13', 'ref.fa'), ('TGT.1', 't1.fa'), ('TGT.2', 't2.fa'),
                                              ('DON.1', 'd1.fa'), ('DON.2', 'd2.fa')])
        errmap = {'TGT.1': '/x/t1.bed', 'DON.2': '/x/d2.bed'}
        options = self._options_with_beds(seqfile, errmap, reference=['CHM13'])
        runs = make_runs(options, GRAPH_EVENT, DOT_HAPLOTYPES)
        self.assertEqual(len(runs), 1)
        run = runs[0]
        # masking is keyed by the seqfile event name (both target and donor), normalized like the seqfile's
        # (TGT.1 is TGT#1); the reference is never masked
        self.assertEqual(run['error_beds'], {'TGT#1': '/x/t1.bed', 'DON#2': '/x/d2.bed'})
        # target report-tagging beds are keyed by graph haplotype (TGT.1 -> hap 1)
        self.assertEqual(run['target_error_beds'], {1: '/x/t1.bed'})

    def test_make_runs_error_beds_reference_free(self):
        seqfile = self._seqfile('rf.seqfile', [('TGT.1', 't1.fa'), ('TGT.2', 't2.fa'), ('DON', 'd.fa')])
        errmap = {'TGT.1': '/x/t1.bed', 'DON': '/x/d.bed'}
        options = self._options_with_beds(seqfile, errmap)
        self.assertTrue(options.referenceFree)
        runs = self._by_name(make_runs(options, GRAPH_EVENT, DOT_HAPLOTYPES))
        # the TGT.1 run (named <base>.hap<N>): the target is renamed to graph_name TGT_1, so its mask
        # bed is keyed by that; the donor keeps its name; the report-tagging bed is keyed by graph hap 0
        r1 = runs['TGT.hap1']
        self.assertEqual(r1['error_beds'], {'DON': '/x/d.bed', 'TGT_1': '/x/t1.bed'})
        self.assertEqual(r1['target_error_beds'], {0: '/x/t1.bed'})
        # the TGT.2 run has no target bed, only the donor mask
        r2 = runs['TGT.hap2']
        self.assertEqual(r2['error_beds'], {'DON': '/x/d.bed'})
        self.assertEqual(r2['target_error_beds'], {})

    def test_target_error_bed_keys(self):
        # the target's BEDs in error_beds (keyed for panpatch) are the ones target_error_beds keys by graph haplotype
        seqfile = self._seqfile('k.seqfile', [('CHM13', 'ref.fa'), ('TGT.1', 't1.fa'), ('TGT.2', 't2.fa'),
                                              ('DON.1', 'd1.fa')])
        errmap = {'TGT.1': '/x/t1.bed', 'TGT.2': '/x/t2.bed', 'DON.1': '/x/d1.bed'}
        run = make_runs(self._options_with_beds(seqfile, errmap, reference=['CHM13']), GRAPH_EVENT, DOT_HAPLOTYPES)[0]
        self.assertEqual(target_error_bed_keys(run), {1: 'TGT#1', 2: 'TGT#2'})
        seqfile = self._seqfile('kf.seqfile', [('TGT.1', 't1.fa'), ('TGT.2', 't2.fa'), ('DON', 'd.fa')])
        errmap = {'TGT.1': '/x/t1.bed', 'DON': '/x/d.bed'}
        runs = self._by_name(make_runs(self._options_with_beds(seqfile, errmap), GRAPH_EVENT, DOT_HAPLOTYPES))
        self.assertEqual(target_error_bed_keys(runs['TGT.hap1']), {0: 'TGT_1'})
        self.assertEqual(target_error_bed_keys(runs['TGT.hap2']), {})       # the donor's BED only
        for run in runs.values():
            self.assertEqual({hap: run['error_beds'][key] for hap, key in target_error_bed_keys(run).items()},
                             run['target_error_beds'])

    # ---- the target's tandem-repeat filter (--maxTargetErrorTandem) ----

    def _fasta(self, name, recs, gz=False, width=60):
        path = os.path.join(self.tempDir, name)
        with (gzip.open(path, 'wt') if gz else open(path, 'w')) as f:
            for rec_name, seq in recs:
                f.write('>{} a description\n'.format(rec_name))
                f.write(''.join(seq[i:i + width] + '\n' for i in range(0, len(seq), width)))
        return path

    def test_error_interval_sequences(self):
        a = ''.join('ACGT'[(i * 7 + i // 3) % 4] for i in range(1000))
        b = 'T' * 50 + 'G' * 50
        for gz in (False, True):
            fa = self._fasta('s.fa.gz' if gz else 's.fa', [('ctgA', a), ('ctgB', b)], gz=gz)
            ivs = [('ctgB', 10, 20), ('ctgA', 100, 200), ('missing', 0, 10), ('ctgA', 990, 1000)]
            # each interval with up to 50 bp either side, clipped to its contig; one on a contig the fasta lacks is skipped
            self.assertEqual(sorted(error_interval_sequences(fa, ivs, context=50)),
                             [(0, 0, b[0:70]), (1, 50, a[50:250]), (3, 940, a[940:1000])])

    def test_covered_fraction(self):
        self.assertEqual(covered_fraction(0, 100, []), 0.0)
        self.assertEqual(covered_fraction(0, 100, [(10, 30), (20, 40), (90, 200)]), 0.4)    # a union, clipped
        self.assertEqual(covered_fraction(50, 100, [(0, 60), (55, 58)]), 0.2)
        self.assertEqual(covered_fraction(5, 5, [(0, 10)]), 0.0)

    def test_tandem_filter_bed(self):
        # a stand-in for FasTAN that calls each run of the unit CATTC a tandem repeat, so the bookkeeping can be checked
        # exactly: record coordinates back through the context, N counted with the repeats, the threshold exclusive
        rng = random.Random(1)
        uniq = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        array = 'CATTC' * 1600
        chr1 = uniq(20000) + array + uniq(20000) + 'N' * 6000 + uniq(26000)    # array [20000, 28000), N [48000, 54000)
        chrJ = uniq(20000) + array + uniq(30000)                                # array [20000, 28000)
        chrK = uniq(20000) + array + uniq(30000)                                # array [20000, 28000)
        chrN = uniq(10000) + 'N' * 4000 + uniq(10000)                           # N [10000, 14000)
        fa = self._fasta('t.fa', [('chr1', chr1), ('chrA', 'A' * 4000 + 'C' * 500), ('chrS', 'ACGT' * 100),
                                  ('chrJ', chrJ), ('chrK', chrK), ('chrN', chrN)])
        bed = self._write('t.bed', '# a comment\n'
                          'chr1\t20000\t28000\n'    # the array: left out
                          'chr1\t30000\t40000\n'    # unique: kept
                          'chr1\t47000\t55000\n'    # 75% N: left out
                          'chrA\t0\t4500\n'         # 89% one base: a tandem repeat of period 1, left out unrun
                          'chrS\t0\t400\n'          # too little sequence to run: judged on its N, kept
                          'missing\t0\t100\n'       # a contig the fasta lacks: kept
                          'chrJ\t20000\t28000\n'    # the array joined to 2 kb beside it: 80% tandem, both left out
                          'chrJ\t28000\t30000\n'
                          'chrJ\t40000\t44000\n'    # unique, joined: kept
                          'chrJ\t44000\t48000\n'
                          'chrK\t28000\t40000\n'    # the array joined to 12 kb beside it (lines out of order): 40%
                          'chrK\t20000\t28000\n'    # tandem, both kept, though the array's own line is all tandem
                          'chrN\t6000\t14000\n')    # 50% N, not more: kept
        records = []

        def fake_fastan(fa_path, work_dir, threads):
            recs = _read_fasta_bytes(fa_path)
            records.extend(recs)
            return {name: [(m.start(), m.end()) for m in re.finditer(rb'(?:CATTC){4,}', seq)] for name, seq in recs.items()}

        out = os.path.join(self.tempDir, 'out.bed')
        n_ivs, dropped = tandem_filter_bed(bed, fa, out, 0.5, self.tempDir, find_tandem=fake_fastan)
        self.assertEqual(n_ivs, 10)
        self.assertEqual([(c, s, e) for c, s, e, _ in dropped], [('chr1', 20000, 28000), ('chr1', 47000, 55000),
                                                                 ('chrA', 0, 4500), ('chrJ', 20000, 30000)])
        self.assertEqual([round(f, 3) for _, _, _, f in dropped], [1.0, 0.75, 1.0, 0.8])
        with open(out) as f:
            self.assertEqual(f.read(), '# a comment\nchr1\t30000\t40000\nchrS\t0\t400\nmissing\t0\t100\n'
                             'chrJ\t40000\t44000\nchrJ\t44000\t48000\nchrK\t28000\t40000\nchrK\t20000\t28000\n'
                             'chrN\t6000\t14000\n')
        self.assertEqual(len(records), 7)       # chrA and chrS were judged without it, and missing not at all
        # a limit of 1 leaves everything in
        self.assertEqual(tandem_filter_bed(bed, fa, out, 1.0, self.tempDir, find_tandem=fake_fastan)[1], [])

    @unittest.skipUnless(shutil.which('FasTAN') and shutil.which('tanbed'), 'FasTAN and tanbed are not on the PATH')
    def test_filter_target_error_beds_fastan(self):
        # FasTAN itself, through the job: an array of a 171 bp unit (as alpha satellite's), its copies ~4% apart, is left
        # out; unique sequence is kept
        rng = random.Random(7)
        uniq = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        unit = uniq(171)
        array = ''.join(''.join(rng.choice('ACGT') if rng.random() < 0.02 else c for c in unit) for _ in range(120))
        fa = self._fasta('t.fa', [('chr1', uniq(30000) + array + uniq(30000))])
        bed = self._write('t.bed', 'chr1\t30000\t{}\nchr1\t5000\t25000\n'.format(30000 + len(array)))
        mode = os.environ.get('CACTUS_BINARIES_MODE')
        os.environ['CACTUS_BINARIES_MODE'] = 'local'
        try:
            out = filter_target_error_beds(self.job, Namespace(maxTargetErrorTandem=0.5), {'name': 'TGT.hap1'},
                                           {0: bed}, {0: fa})
        finally:
            if mode is None:
                del os.environ['CACTUS_BINARIES_MODE']
            else:
                os.environ['CACTUS_BINARIES_MODE'] = mode
        with open(out[0]) as f:
            self.assertEqual(f.read(), 'chr1\t5000\t25000\n')

    def test_make_runs_rejects_unmatched_event(self):
        seqfile = self._seqfile('u.seqfile', [('CHM13', 'ref.fa'), ('TGT.1', 't1.fa'), ('DON.1', 'd1.fa')])
        errmap = {'BOGUS': '/x/bogus.bed'}
        options = self._options_with_beds(seqfile, errmap, reference=['CHM13'])
        with self.assertRaises(RuntimeError):
            make_runs(options, GRAPH_EVENT, DOT_HAPLOTYPES)

    def test_make_runs_rejects_unmatched_event_sharing_a_bed_path(self):
        # the typo check must compare EVENT names, not bed paths: a mistyped event that happens to share
        # a bed file with a valid event must still be caught (else the mistyped assembly is silently
        # never masked -- and if it was a donor, its error sequence can leak into a patch)
        seqfile = self._seqfile('sp.seqfile', [('CHM13', 'ref.fa'), ('TGT.1', 't1.fa'), ('DON.1', 'd1.fa')])
        shared = '/x/shared.bed'
        errmap = {'TGT.1': shared, 'BOGUS': shared}
        options = self._options_with_beds(seqfile, errmap, reference=['CHM13'])
        with self.assertRaises(RuntimeError):
            make_runs(options, GRAPH_EVENT, DOT_HAPLOTYPES)
        # naming the reference (which is never masked) is likewise an error, even sharing a valid path
        errmap2 = {'DON.1': shared, 'CHM13': shared}
        options2 = self._options_with_beds(seqfile, errmap2, reference=['CHM13'])
        with self.assertRaises(RuntimeError):
            make_runs(options2, GRAPH_EVENT, DOT_HAPLOTYPES)


    # ---- batch summary (summarize_report / write_batch_summary) ----

    def _report(self, name, rows, contigs, summary=()):
        """ write a per-sample panpatch report (header + gap_origin column + rows + #Summary lines + #Contig
        cap lines).  rows: (type, decision, gap_origin) tuples; contigs: (left_yes, right_yes) tuples """
        hdr = _REPORT_HEADER + '\tgap_origin'
        def row(typ, decision, gap_origin):
            return '\t'.join(['chr1', '1', typ, 'TGT#1#chr1', '100', 'DON#1#d', '20', '0',
                              '.', '100.0', '100.0', decision, '.', '10', '20', gap_origin])
        lines = [hdr] + [row(*r) for r in rows] + list(summary)
        for i, (l, r) in enumerate(contigs):
            lines.append('#Contig TGT#1#c{} len=90000bp left={}(0.9) right={}(0.9)'.format(
                i, 'YES' if l else 'NO', 'YES' if r else 'NO'))
        return self._write(name, '\n'.join(lines) + '\n')

    def test_summarize_report(self):
        rep = self._report('a.tsv',
                           [('gap-fill', 'accepted', 'bed-gap'),
                            ('gap-fill', 'accepted', 'N-gap'),
                            ('scaffold', 'accepted', '.'),
                            ('telomere', 'accepted', '.'),
                            ('telomere', 'accepted', 'bed-gap'),   # a telomere patch over a flagged tip counts too
                            ('gap-fill', 'rejected', 'N-gap'),     # rejected -> not counted
                            ('passthrough', 'accepted', '.')],     # passthrough -> not a patch
                           [(True, True), (True, False), (False, True)])  # only the first is T2T
        s = summarize_report(rep)
        self.assertEqual(s['counts'], {'gap-fill': 2, 'bed-gap': 2, 'scaffold': 1, 'telomere': 2})
        self.assertEqual(s['t2t'], 1)
        self.assertEqual(s['contigs'], 3)
        self.assertEqual(s['rejected'], 1)               # the rejected gap-fill; passthrough is neither
        self.assertEqual(s['output_records'], 3)         # no #Summary line: counted from the cap lines

    def test_t2t_needs_no_gap(self):
        # telomere-to-telomere = a telomere at both ends and no N gap; a cap line from a panpatch too old to
        # print gaps= is judged by its telomeres alone
        line = '#Contig TGT#1#c len=90000bp left={}(0.9) right={}(0.9){}'
        self.assertTrue(contig_is_t2t(line.format('YES', 'YES', ' gaps=0')))
        self.assertFalse(contig_is_t2t(line.format('YES', 'YES', ' gaps=1')))
        self.assertFalse(contig_is_t2t(line.format('YES', 'YES', ' gaps=12')))
        self.assertFalse(contig_is_t2t(line.format('YES', 'NO', ' gaps=0')))
        self.assertTrue(contig_is_t2t(line.format('YES', 'YES', '')))
        rep = self._write('g.tsv', _REPORT_HEADER + '\tgap_origin\n'
                          + '\n'.join(line.format('YES', 'YES', g) for g in (' gaps=0', ' gaps=1', ' gaps=0')) + '\n')
        self.assertEqual(summarize_report(rep)['t2t'], 2)

    def test_summarize_report_summary_lines(self):
        # gather_panpatch's #Summary lines (one per haplotype) give the before/after counts, summed
        rep = self._report('s.tsv', [('telomere', 'accepted', '.')], [(True, True)],
                           summary=['#Summary hap=1 input_contigs=30 input_t2t=19 output_records=30 t2t=22',
                                    '#Summary hap=2 input_contigs=24 input_t2t=19 output_records=25 t2t=23'])
        s = summarize_report(rep)
        self.assertEqual((s['input_contigs'], s['input_t2t'], s['output_records'], s['t2t']), (54, 38, 55, 45))

    def test_concat_reports_replaces_stale_cap_lines(self):
        # two per-chromosome reports, each carrying panpatch's own per-path cap lines (which describe
        # contributing graph paths before the target revert, donors included), plus the cap lines the
        # caller measured on the finished assembly.  The concatenation must keep one header, keep every
        # patch row, drop ALL of panpatch's cap lines, and end with exactly the measured ones -- so a
        # summary built from it counts output contigs and post-revert telomeres.
        hdr = _REPORT_HEADER + '\tgap_origin'
        def row(chrom, typ):
            return '\t'.join([chrom, '1', typ, 'TGT#1#c', '100', 'DON#1#d', '20', '0',
                              '.', '100.0', '100.0', 'accepted', '.', '10', '20', 'N-gap'])
        c1 = self._write('chr1.tsv', hdr + '\n' + row('chr1', 'gap-fill') + '\n'
                         + '#Contig TGT#1#c#0 len=100bp left=NO(0.001) right=YES(0.9)\n'
                         + '#Contig DON#1#d#0 len=20bp left=NO(0.0) right=NO(0.0)\n')
        c2 = self._write('chr2.tsv', hdr + '\n' + row('chr2', 'telomere') + '\n'
                         + '#Contig TGT#1#c2#0 len=50bp left=NO(0.0) right=NO(0.0)\n')
        measured = ['#Contig chr1_hap_1 len=120bp left=YES(0.95) right=YES(0.93)',
                    '#Contig chr2_hap_1 len=60bp left=YES(0.91) right=NO(0.02)']
        out = os.path.join(self.tempDir, 'joined.tsv')
        concat_reports([c1, c2], measured, out)

        lines = [l.rstrip('\n') for l in open(out)]
        self.assertEqual(lines.count(hdr), 1)                       # header exactly once
        self.assertIn(row('chr1', 'gap-fill'), lines)                # every patch row survives
        self.assertIn(row('chr2', 'telomere'), lines)
        caps = [l for l in lines if l.startswith('#Contig')]
        self.assertEqual(caps, measured)                             # only the measured caps remain
        self.assertNotIn('DON#1#d#0', '\n'.join(caps))               # donor path line is gone

        # and the summary built from it reports output contigs / post-revert T2T
        s = summarize_report(out)
        self.assertEqual(s['contigs'], 2)
        self.assertEqual(s['t2t'], 1)

    def test_write_batch_summary(self):
        a = self._report('A.tsv',
                         [('gap-fill', 'accepted', 'bed-gap'), ('gap-fill', 'accepted', 'N-gap'),
                          ('scaffold', 'accepted', '.')],
                         [(True, True), (True, True), (False, True)],          # 2 T2T of 3
                         summary=['#Summary hap=1 input_contigs=4 input_t2t=1 output_records=3 t2t=2'])
        b = self._report('B.tsv',
                         [('gap-fill', 'accepted', 'N-gap'), ('telomere', 'accepted', '.'),
                          ('gap-fill', 'rejected', 'N-gap')],
                         [(True, True), (False, False)])                       # 1 T2T of 2, no #Summary
        summaries = [('A', {'A.tsv': a}), ('B', {'B.tsv': b}),
                     ('SKIP', None), ('SKIP2', {'other.fa': a})]               # runs with no report -> skipped
        options = Namespace(outDir=self.tempDir)
        write_batch_summary(self.job, options, summaries)

        with open(os.path.join(self.tempDir, 'panpatch-summary.tsv')) as f:
            lines = [l.rstrip('\n') for l in f if not l.startswith('#')]
        table = {c.split('\t')[0]: c.split('\t')[1:] for c in lines[1:]}   # skip the column header
        self.assertEqual(lines[0], 'run\tinput_contigs\toutput_records\tinput_t2t\tt2t\tpatches\tgap_fill\tbed_gap\t'
                                   'scaffold\ttelomere\trejected')
        self.assertEqual(table['A'], ['4', '3', '1', '2', '3', '2', '1', '1', '0', '0'])
        self.assertEqual(table['B'], ['0', '2', '0', '1', '2', '1', '0', '0', '1', '1'])
        self.assertEqual(table['TOTAL'], ['4', '5', '1', '3', '5', '3', '1', '1', '1', '1'])
        self.assertNotIn('SKIP', table)                                    # runs without a report are omitted

    # ---- the delivered names, and the completeness check ----

    def test_clean_path_name(self):
        renames = run_renames({'referenceFree': True, 'samples': ['HG2_v_1', 'DON'], 'target_sample': 'HG2.v', 'hap': '1'})
        self.assertEqual(renames, {'HG2_v_1': ('HG2.v', '1')})
        self.assertEqual(clean_path_name('HG2_v_1#0#chr8_haplotype1-0000014', renames), 'HG2.v#1#chr8_haplotype1-0000014')
        self.assertEqual(clean_path_name('DON#2#CM1.1#0', renames), 'DON#2#CM1.1')           # phase block dropped
        self.assertEqual(clean_path_name('TGT#1#c#0[8-100]', {}), 'TGT#1#c[8-100]')         # a subrange is kept
        self.assertEqual(clean_path_name('chr8_hap_1', renames), 'chr8_hap_1')              # not a path: as it is
        self.assertEqual(run_renames({'referenceFree': False}), {})

    def test_clean_report_and_bed_names(self):
        renames = {'T_1': ('T', '1')}
        rep = self._write('r.tsv', _REPORT_HEADER + '\tgap_origin\n' + '\t'.join(
            ['c1', '0', 'gap-fill', 'T_1#0#c1', '100', 'DON#2#d#0', '20', '5', '.', '.', '.', 'rejected',
             'contig reverted: contig T_1#0#c1 scaffolded to DON#2#d#0', '10', '15', 'N-gap']) + '\n'
            + '#Contig c1 len=100bp left=YES(0.9) right=YES(0.9)\n')
        clean_report_names(rep, renames, hap='1')
        f = open(rep).read().split('\n')[1].split('\t')
        self.assertEqual((f[1], f[3], f[5]), ('1', 'T#1#c1', 'DON#2#d'))
        self.assertEqual(f[12], 'contig reverted: contig T#1#c1 scaffolded to DON#2#d')
        self.assertIn('#Contig c1 len=100bp', open(rep).read())                         # cap lines untouched
        bed = self._write('b.bed', '#Patched assembly on c1 for T_1#0:\nT_1#0#c1\t0\t10\t+\nDON#2#d#0\t3\t9\t+\n\n'
                                   '#Unpatched (input contig) on c2 for T_1#0:\nT_1#0#c2\t0\t50\t+\n')
        clean_bed_names(bed, renames)
        self.assertEqual(open(bed).read(), '#Patched assembly on c1 for T#1:\nT#1#c1\t0\t10\t+\nDON#2#d\t3\t9\t+\n\n'
                                           '#Unpatched (input contig) on c2 for T#1:\nT#1#c2\t0\t50\t+\n')

    def test_output_record_names(self):
        # reference-free: every record is one input contig and takes its name
        self.assertEqual(output_record_name('c1_hap_0', True), 'c1')
        self.assertEqual(output_record_name('T_1#0#c2', True), 'c2')
        # reference-based: a patched chromosome keeps <chrom>_hap_<N>; a contig written as it is, its name
        self.assertEqual(output_record_name('chr8_hap_1', False), 'chr8_hap_1')
        self.assertEqual(output_record_name('T#1#c3#0', False), 'c3')
        fa = self._write('n.fa', '>c1_hap_0\nACGT\n>T_1#0#c2\nGG\n')
        self.assertEqual(rename_fasta_records(fa, True), ['c1', 'c2'])
        self.assertEqual(open(fa).read(), '>c1\nACGT\n>c2\nGG\n')
        dup = self._write('d.fa', '>c1_hap_0\nA\n>T_1#0#c1\nC\n')                        # both would be c1
        with self.assertRaises(RuntimeError):
            rename_fasta_records(dup, True)

    def test_missing_target_contigs(self):
        bed = self._write('m.bed', '#Patched assembly on chr1 for T#1:\nT#1#a#0\t0\t10\t+\nDON#1#x#0\t0\t5\t+\n'
                                   'T#1#b#0\t0\t9\t+\n\n'
                                   '#Patched assembly on chr1 for T#2:\nT#2#z#0\t0\t10\t+\n\n')
        rep = self._write('m.tsv', _REPORT_HEADER + '\n' + '\t'.join(
            ['chr1', '1', 'telomere', 'T#1#c#0', '9', 'DON#1#x#0', '20', '9', '.', '.', '.', 'accepted', '.', '0', '9']) + '\n')
        # a and b are in a patched record, c was taken over by an accepted patch, d is written as itself,
        # e is nowhere (and z is the other haplotype's, so it does not count for hap 1)
        missing = missing_target_contigs(['a', 'b', 'c', 'd', 'e', 'z'], ['chr1_hap_1', 'T#1#d#0'], bed, rep, 'T', 1)
        self.assertEqual(missing, ['e', 'z'])
    def test_rescue_dropped_contigs(self):
        # contigs no graph saw are appended as they are, as SAMPLE#HAP#CONTIG, and returned with their lengths
        target = self._write('t.fa', '>c1 a comment\nACGT\nAC\n>c2\nGG\n>c3\nTTTTT\n')
        out = self._write('o.fa', '>c1_hap_0\nACGTAC\n')
        rescued = rescue_dropped_contigs(self.job, out, target, {'0#c1'}, 'T_1', 0, label='T.hap1')
        self.assertEqual(rescued, [('c2', 2), ('c3', 5)])
        self.assertEqual(open(out).read(), '>c1_hap_0\nACGTAC\n>T_1#0#c2\nGG\n>T_1#0#c3\nTTTTT\n')

    def test_report_columns_match(self):
        # the header gather writes for a report panpatch never wrote is panpatch's own
        self.assertEqual(PANPATCH_REPORT_COLUMNS, _REPORT_HEADER)

if __name__ == '__main__':
    unittest.main()
