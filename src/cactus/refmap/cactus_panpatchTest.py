#!/usr/bin/env python3

"""
Unit tests for the pure (non-Toil) logic of cactus-panpatch: the fasta-header sanitizer that has to
match cactus_sanitizeFastaHeaders, the sample-name helpers, and the run resolution that turns a
seqfile (or chromfile) plus options into the list of cactus-pangenome+panpatch runs.

These are fast and offline: make_runs only reads the seqfile text (SeqFile does not open the fastas),
so the fasta paths below are dummies.  The end-to-end pipeline is covered by evolverTest.py.
"""

import os
import gzip
import unittest
import tempfile
import shutil
from argparse import Namespace

from cactus.refmap.cactus_panpatch import sanitize_contig_name, sample_base, sample_hap
from cactus.refmap.cactus_panpatch import make_runs, panpatch_validate_options
from cactus.refmap.cactus_panpatch import read_error_bed_manifest, parse_error_bed, intervals_overlap, revcomp
from cactus.refmap.cactus_panpatch import parse_bed_blocks, mask_assembly_errors, revert_target_to_original
from cactus.refmap.cactus_panpatch import tag_report_gap_origin

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

    def test_revcomp(self):
        self.assertEqual(revcomp(b'ACGTN'), b'NACGT')
        self.assertEqual(revcomp(b'acgtn'), b'nacgt')          # case preserved
        self.assertEqual(revcomp(b'AAAC'), b'GTTT')
        self.assertEqual(revcomp(revcomp(b'ACGTRYSWKM')), b'ACGTRYSWKM')  # IUPAC round-trips

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

    # ---- mask_assembly_errors ----

    def test_mask_assembly_errors_plain(self):
        t = 'A' * 30 + 'C' * 20 + 'G' * 50            # 100 bp, one contig
        fa = self._write('in.fa', '>chr1 some description\n{}\n>chr2\n{}\n'.format(t, 'T' * 40))
        bed = self._write('m.bed', 'chr1\t30\t50\n')  # matches the header's FIRST token only
        out = mask_assembly_errors(self.job, fa, bed)
        recs = _read_fasta_bytes(out)
        # mask preserves the header verbatim (the description survives; sanitize truncates it later), but
        # matches the bed against the header's FIRST token only
        self.assertIn('chr1 some description', recs)
        masked = recs['chr1 some description']
        self.assertEqual(len(masked), 100)                              # length preserved
        self.assertEqual(masked, b'A' * 30 + b'N' * 20 + b'G' * 50)
        self.assertEqual(recs['chr2'], b'T' * 40)                       # untouched contig unchanged

    def test_mask_assembly_errors_gzip_input(self):
        t = 'A' * 100
        raw = '>chr1\n{}\n'.format(t).encode()
        gz = os.path.join(self.tempDir, 'in.fa.gz')
        with gzip.open(gz, 'wb') as f:
            f.write(raw)
        bed = self._write('m.bed', 'chr1\t10\t20\n')
        out = mask_assembly_errors(self.job, gz, bed)
        # gzipped input is read transparently; output is plaintext (the sanitize step downstream sniffs
        # magic bytes, so there is no reason to recompress)
        with open(out, 'rb') as f:
            self.assertNotEqual(f.read(2), b'\x1f\x8b')
        recs = _read_fasta_bytes(out)
        self.assertEqual(recs['chr1'], b'A' * 10 + b'N' * 10 + b'A' * 80)

    # ---- revert_target_to_original ----

    def test_revert_restores_target_keeps_donor(self):
        # pristine target haplotype (what the revert reaches back to)
        T   = b'A' * 30 + b'C' * 20 + b'G' * 10 + b'T' * 20 + b'A' * 20   # chr1, 100 bp
        T2  = b'ACGTAC' * 4                                               # chr2, 24 bp
        T1b = b'TTTTGGGGCCCCAAAA'                                         # chr1b, 16 bp
        pristine = self._write('pristine.fa',
                               '>chr1\n{}\n>chr2\n{}\n>chr1b\n{}\n'.format(T.decode(), T2.decode(), T1b.decode()))

        # the panpatch output for hap 1: a masked-but-unpatched gap survives as N inside a patched
        # record (chr1: target 0-60 carries N at 30-50, donor fills 60-80, target 80-100); chr2 is a
        # reverse-strand patched record with a masked gap; chr1b is a reverted contig full of masked Ns;
        # a donor contig record must pass through untouched
        chr1_out = (T[0:30] + b'N' * 20 + T[50:60]) + (b'g' * 20) + T[80:100]
        chr2_masked = T2[0:4] + b'N' * 4 + T2[8:24]
        chr2_out = revcomp(chr2_masked)
        chr1b_out = T1b[0:4] + b'N' * 4 + T1b[8:16]
        don_out = b'g' * 12
        combined = self._write('hap1.fa',
                               '>chr1_hap_1\n{}\n>chr2_hap_1\n{}\n>TGT#1#chr1b\n{}\n>DON#1#dcontig2\n{}\n'.format(
                                   chr1_out.decode(), chr2_out.decode(), chr1b_out.decode(), don_out.decode()))

        bed = ('#Patched assembly on chr1 for TGT#1:\n'
               'TGT#1#chr1\t0\t60\t+\n'
               'DON#1#dcontig\t60\t80\t+\n'
               'TGT#1#chr1\t80\t100\t+\n'
               '\n'
               '#Patched assembly on chr2 for TGT#1:\n'
               'TGT#1#chr2\t0\t24\t-\n')
        blocks = parse_bed_blocks(self._write('hap1.bed', bed))

        revert_target_to_original(self.job, combined, blocks, pristine, 'TGT', 1)
        recs = _read_fasta_bytes(combined)

        # chr1: the N gap is restored from pristine, the donor bytes are kept verbatim
        self.assertEqual(recs['chr1_hap_1'], T[0:60] + b'g' * 20 + T[80:100])
        self.assertEqual(recs['chr1_hap_1'][30:50], b'C' * 20)   # was N, now original
        self.assertEqual(recs['chr1_hap_1'][60:80], b'g' * 20)   # donor untouched
        # chr2: reverse-strand target fully restored (revcomp of pristine), no Ns
        self.assertEqual(recs['chr2_hap_1'], revcomp(T2))
        self.assertNotIn(b'N', recs['chr2_hap_1'])
        # chr1b: whole reverted contig restored from pristine
        self.assertEqual(recs['TGT#1#chr1b'], T1b)
        # donor contig record passes through unchanged
        self.assertEqual(recs['DON#1#dcontig2'], don_out)

    def test_revert_cursor_mismatch_raises(self):
        # a bed block whose interval lengths do not sum to the record length is a corruption we must
        # not silently accept (it would mean the byte-exact contract was violated)
        pristine = self._write('p.fa', '>chr1\n{}\n'.format('A' * 100))
        combined = self._write('h.fa', '>chr1_hap_1\n{}\n'.format('A' * 50))   # record is 50 bp
        bed = '#Patched assembly on chr1 for TGT#1:\nTGT#1#chr1\t0\t100\t+\n'  # claims 100 bp
        blocks = parse_bed_blocks(self._write('h.bed', bed))
        with self.assertRaises(RuntimeError):
            revert_target_to_original(self.job, combined, blocks, pristine, 'TGT', 1)

    # ---- tag_report_gap_origin ----

    def test_tag_report_gap_origin(self):
        # a target error bed on hap 1, chr1, covering [30,50)
        err = self._write('t1.bed', 'chr1\t30\t50\n')
        target_error_bed_ids = {1: err}

        def row(chrom, hap, typ, target, ts, te):
            return '\t'.join([chrom, str(hap), typ, target, '100', 'DON#1#d', '20', '0',
                              '.', '100.0', '100.0', 'accepted', '.', str(ts), str(te)])
        report = self._write('report.tsv', _REPORT_HEADER + '\n'
                             + row('chr1', 1, 'gap-fill', 'TGT#1#chr1', 35, 45) + '\n'    # overlaps -> bed-gap
                             + row('chr1', 1, 'gap-fill', 'TGT#1#chr1', 70, 80) + '\n'    # no overlap -> N-gap
                             + row('chr1', 1, 'telomere', 'TGT#1#chr1', 0, 10) + '\n'     # not gap-fill -> .
                             + row('chr1', 2, 'gap-fill', 'TGT#2#chr1', 35, 45) + '\n'    # hap 2 has no bed -> N-gap
                             + '#Contig TGT#1#chr1 len=100\n')

        tag_report_gap_origin(self.job, report, target_error_bed_ids)
        with open(report) as f:
            lines = [l.rstrip('\n') for l in f]

        self.assertEqual(lines[0], _REPORT_HEADER + '\tgap_origin')
        self.assertTrue(lines[1].endswith('\tbed-gap'))
        self.assertTrue(lines[2].endswith('\tN-gap'))
        self.assertTrue(lines[3].endswith('\t.'))          # telomere row
        self.assertTrue(lines[4].endswith('\tN-gap'))      # hap 2, no error bed
        self.assertEqual(lines[5], '#Contig TGT#1#chr1 len=100')   # comment passes through, no column

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
        runs = make_runs(options, GRAPH_EVENT)
        self.assertEqual(len(runs), 1)
        run = runs[0]
        # masking is keyed by the seqfile event name (both target and donor); the reference is never masked
        self.assertEqual(run['mask_beds'], {'TGT.1': '/x/t1.bed', 'DON.2': '/x/d2.bed'})
        # target report-tagging beds are keyed by graph haplotype (TGT.1 -> hap 1)
        self.assertEqual(run['target_error_beds'], {1: '/x/t1.bed'})

    def test_make_runs_error_beds_reference_free(self):
        seqfile = self._seqfile('rf.seqfile', [('TGT.1', 't1.fa'), ('TGT.2', 't2.fa'), ('DON', 'd.fa')])
        errmap = {'TGT.1': '/x/t1.bed', 'DON': '/x/d.bed'}
        options = self._options_with_beds(seqfile, errmap)
        self.assertTrue(options.referenceFree)
        runs = self._by_name(make_runs(options, GRAPH_EVENT))
        # the TGT.1 run (named <base>.hap<N>): the target is renamed to graph_name TGT_1, so its mask
        # bed is keyed by that; the donor keeps its name; the report-tagging bed is keyed by graph hap 0
        r1 = runs['TGT.hap1']
        self.assertEqual(r1['mask_beds'], {'DON': '/x/d.bed', 'TGT_1': '/x/t1.bed'})
        self.assertEqual(r1['target_error_beds'], {0: '/x/t1.bed'})
        # the TGT.2 run has no target bed, only the donor mask
        r2 = runs['TGT.hap2']
        self.assertEqual(r2['mask_beds'], {'DON': '/x/d.bed'})
        self.assertEqual(r2['target_error_beds'], {})

    def test_make_runs_rejects_unmatched_event(self):
        seqfile = self._seqfile('u.seqfile', [('CHM13', 'ref.fa'), ('TGT.1', 't1.fa'), ('DON.1', 'd1.fa')])
        errmap = {'BOGUS': '/x/bogus.bed'}
        options = self._options_with_beds(seqfile, errmap, reference=['CHM13'])
        with self.assertRaises(RuntimeError):
            make_runs(options, GRAPH_EVENT)

    def test_make_runs_rejects_unmatched_event_sharing_a_bed_path(self):
        # the typo check must compare EVENT names, not bed paths: a mistyped event that happens to share
        # a bed file with a valid event must still be caught (else the mistyped assembly is silently
        # never masked -- and if it was a donor, its error sequence can leak into a patch)
        seqfile = self._seqfile('sp.seqfile', [('CHM13', 'ref.fa'), ('TGT.1', 't1.fa'), ('DON.1', 'd1.fa')])
        shared = '/x/shared.bed'
        errmap = {'TGT.1': shared, 'BOGUS': shared}
        options = self._options_with_beds(seqfile, errmap, reference=['CHM13'])
        with self.assertRaises(RuntimeError):
            make_runs(options, GRAPH_EVENT)
        # naming the reference (which is never masked) is likewise an error, even sharing a valid path
        errmap2 = {'DON.1': shared, 'CHM13': shared}
        options2 = self._options_with_beds(seqfile, errmap2, reference=['CHM13'])
        with self.assertRaises(RuntimeError):
            make_runs(options2, GRAPH_EVENT)


if __name__ == '__main__':
    unittest.main()
