#!/usr/bin/env python3

"""
Unit tests for the pure (non-Toil) logic behind the pangenome exclusion report: the hal2vg name
resolution that has to match hal2vg.cpp, PanSN path parsing, interval algebra, and the parsers for
the baseline table, the split log and vg depth output.

These are fast and offline.  The end-to-end report is covered by evolverTest.py.
"""

import gzip
import subprocess
import os
import shutil
import tempfile
import unittest

from cactus.refmap.pangenome_exclusions import (
    compute_exclusions, sublabel_dropped, write_summary, merge_sorted_bed_stream, pansn_prefix_to_event,
    DEPTH_TO_BED_AWK, PHASE_REASON, check_baseline_header, detect_problems,
    BASELINE_HEADER, baseline_by_key, contig_sizes_from_fai,
    event_by_pansn_prefix, event_to_pansn_prefix, intersect_intervals, merge_intervals,
    parse_path_name, parse_split_log, read_baseline_tsv, resolve_subpath_naming,
    safe_event_filename, subtract_intervals, total_bp, write_baseline_tsv)


class TestSubpathNaming(unittest.TestCase):

    def test_plain_contig(self):
        self.assertEqual(resolve_subpath_naming('chr1'), ('chr1', 0, None))

    def test_single_subpath(self):
        self.assertEqual(resolve_subpath_naming('chr1_sub_999_2000'), ('chr1', 999, 1001))

    def test_nested_subpath(self):
        # offsets sum; the length is the innermost (first-parsed) one
        self.assertEqual(resolve_subpath_naming('a_sub_10_20_sub_2_5'), ('a', 12, 3))

    def test_non_numeric_suffix_is_opaque(self):
        # a contig genuinely called "scaf_sub_unit" is not a coordinate encoding
        self.assertEqual(resolve_subpath_naming('scaf_sub_unit'), ('scaf_sub_unit', 0, None))

    def test_zero_offset_subpath(self):
        self.assertEqual(resolve_subpath_naming('chrI_sub_0_219929'), ('chrI', 0, 219929))


class TestEventToPansn(unittest.TestCase):

    def test_haploid(self):
        self.assertEqual(event_to_pansn_prefix('S288C'), 'S288C#0')

    def test_diploid(self):
        self.assertEqual(event_to_pansn_prefix('HG002.1'), 'HG002#1')
        self.assertEqual(event_to_pansn_prefix('HG002.2'), 'HG002#2')

    def test_non_numeric_suffix_is_not_a_haplotype(self):
        self.assertEqual(event_to_pansn_prefix('HG002.verkko'), 'HG002.verkko#0')

    def test_numeric_suffix_after_dotted_name(self):
        self.assertEqual(event_to_pansn_prefix('HG002.verkko.1'), 'HG002.verkko#1')


class TestParsePathName(unittest.TestCase):

    def test_reference_sense(self):
        self.assertEqual(parse_path_name('S288C#0#chrI', 219929), ('S288C#0', 'chrI', 0, 219929))

    def test_haplotype_sense(self):
        self.assertEqual(parse_path_name('Y12#0#chrI#0', 197190), ('Y12#0', 'chrI', 0, 197190))

    def test_subrange(self):
        self.assertEqual(parse_path_name('SK1#0#chrI#0[14059-228861]', 214802),
                         ('SK1#0', 'chrI', 14059, 228861))

    def test_subrange_width_equals_path_length(self):
        # the invariant the coverage job asserts on every line
        _p, _c, start, end = parse_path_name('SK1#0#chrI#0[14059-228861]', 214802)
        self.assertEqual(end - start, 214802)

    def test_minigraph_path_still_parses(self):
        self.assertEqual(parse_path_name('_MINIGRAPH_#0#s113[0-405]', 405),
                         ('_MINIGRAPH_#0', 's113', 0, 405))

    def test_non_pansn_is_orphan(self):
        self.assertIsNone(parse_path_name('chrI', 100))
        self.assertIsNone(parse_path_name('sample#0', 100))


class TestIntervals(unittest.TestCase):

    def test_merge(self):
        self.assertEqual(merge_intervals([(5, 10), (0, 3), (2, 6)]), [(0, 10)])
        self.assertEqual(merge_intervals([(0, 3), (3, 6)]), [(0, 6)])
        self.assertEqual(merge_intervals([(0, 3), (4, 6)]), [(0, 3), (4, 6)])
        self.assertEqual(merge_intervals([]), [])
        self.assertEqual(merge_intervals([(5, 5)]), [])

    def test_subtract(self):
        self.assertEqual(subtract_intervals([(0, 100)], [(10, 20)]), [(0, 10), (20, 100)])
        self.assertEqual(subtract_intervals([(0, 100)], [(0, 100)]), [])
        self.assertEqual(subtract_intervals([(0, 100)], []), [(0, 100)])
        self.assertEqual(subtract_intervals([], [(0, 10)]), [])
        self.assertEqual(subtract_intervals([(0, 10), (20, 30)], [(5, 25)]), [(0, 5), (25, 30)])
        self.assertEqual(subtract_intervals([(0, 100)], [(0, 10), (90, 100)]), [(10, 90)])

    def test_intersect(self):
        self.assertEqual(intersect_intervals([(0, 100)], [(10, 20)]), [(10, 20)])
        self.assertEqual(intersect_intervals([(0, 10)], [(20, 30)]), [])
        self.assertEqual(intersect_intervals([(0, 10), (20, 30)], [(5, 25)]), [(5, 10), (20, 25)])

    def test_total_bp(self):
        self.assertEqual(total_bp([(0, 10), (20, 25)]), 15)

    def test_phase_subtraction_is_nested(self):
        # the core of the report: missing_clip - missing_full is what the clip phase removed
        baseline = [(0, 1000)]
        cov_full = [(0, 1000)]
        cov_clip = [(100, 900)]
        missing_full = subtract_intervals(baseline, cov_full)
        missing_clip = subtract_intervals(baseline, cov_clip)
        self.assertEqual(missing_full, [])
        self.assertEqual(subtract_intervals(missing_clip, missing_full),
                         [(0, 100), (900, 1000)])


class TestHprcTwoFragmentClosure(unittest.TestCase):
    """ an input contig split into two _sub_ fragments must close to zero missing.  without the
    hal2vg name resolution this reports the whole contig lost AND the whole contig as orphan
    coverage, which is the single most damaging way this feature could be wrong """

    def test_two_fragments_close_to_zero(self):
        rows = [('HG02080.1', 'HG02080#1', 'JAHEOW010000073.1', 0, 7238466),
                ('HG02080.1', 'HG02080#1', 'JAHEOW010000073.1', 7238466, 5630658)]
        baseline = baseline_by_key(rows)
        key = ('HG02080#1', 'JAHEOW010000073.1')
        self.assertEqual(baseline[key], [(0, 12869124)])

        coverage = []
        for name, length in [('HG02080#1#JAHEOW010000073.1#0[0-7238466]', 7238466),
                             ('HG02080#1#JAHEOW010000073.1#0[7238466-12869124]', 5630658)]:
            pansn, contig, start, end = parse_path_name(name, length)
            self.assertEqual((pansn, contig), key)
            self.assertEqual(end - start, length)
            coverage.append((start, end))

        self.assertEqual(subtract_intervals(baseline[key], coverage), [])

    def test_fai_round_trip(self):
        tmp_dir = tempfile.mkdtemp()
        try:
            fai_path = os.path.join(tmp_dir, 'x.fa.fai')
            with open(fai_path, 'w') as fai_file:
                fai_file.write('id=HG02080.1|JAHEOW010000073.1_sub_0_7238466\t7238466\t0\t60\t61\n')
                fai_file.write('id=HG02080.1|JAHEOW010000073.1_sub_7238466_12869124\t5630658\t0\t60\t61\n')
            rows = contig_sizes_from_fai(fai_path, 'HG02080.1')
            self.assertEqual(len(rows), 2)
            self.assertEqual(baseline_by_key(rows)[('HG02080#1', 'JAHEOW010000073.1')],
                             [(0, 12869124)])
        finally:
            shutil.rmtree(tmp_dir)

    def test_fai_strips_foreign_id_prefix(self):
        # the sanitizer passes through a pre-existing id=X| even when X is not the event
        tmp_dir = tempfile.mkdtemp()
        try:
            fai_path = os.path.join(tmp_dir, 'x.fa.fai')
            with open(fai_path, 'w') as fai_file:
                fai_file.write('id=SOMETHINGELSE|chrI\t219929\t0\t60\t61\n')
            rows = contig_sizes_from_fai(fai_path, 'S288C')
            self.assertEqual(rows[0][1], 'S288C#0')
            self.assertEqual(rows[0][2], 'chrI')
        finally:
            shutil.rmtree(tmp_dir)


class TestBaselineTsv(unittest.TestCase):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def test_round_trip(self):
        rows = [('S288C', 'S288C#0', 'chrI', 0, 219929)]
        path = os.path.join(self.tmp_dir, 'sizes.tsv')
        write_baseline_tsv(rows, path)
        self.assertEqual(read_baseline_tsv(path), rows)

    def test_rejects_contig_sizes_matrix(self):
        """ chrom-subproblems/contig_sizes.tsv lives in the same output tree and has 10 columns of
        completely different meaning.  parsing it would yield a silent 100%-loss report """
        path = os.path.join(self.tmp_dir, 'contig_sizes.tsv')
        with open(path, 'w') as out_file:
            out_file.write('Contig\tDBVPG6044\tS288C\tSK1\tUWOPS034614\tY12\tYPS128\t'
                           '_MINIGRAPH_\tmin\tmax\tavg\n')
            out_file.write('chrXV\t1062383\t1091343\t1053869\t1062502\t1048295\t1072763\t'
                           '1091343\t1048295\t1091343\t1068928\n')
        with self.assertRaises(RuntimeError) as caught:
            read_baseline_tsv(path)
        message = str(caught.exception)
        self.assertIn('contig_sizes.tsv', message)
        self.assertIn('input-contig-sizes.tsv', message)

    def test_event_lookup(self):
        rows = [('S288C', 'S288C#0', 'chrI', 0, 219929),
                ('HG002.1', 'HG002#1', 'chr1', 0, 100)]
        self.assertEqual(event_by_pansn_prefix(rows), {'S288C#0': 'S288C', 'HG002#1': 'HG002.1'})


class TestDepthToBed(unittest.TestCase):
    """ the vg depth -> BED conversion, which is a one-line awk program in the pipe (so that a whole
    genome's per-base output never reaches python) plus a streaming merge.  the awk program is
    pinned here by running it for real: an off-by-one would shift every reference gap """

    @staticmethod
    def depth_lines(path, depths):
        """ vg depth prints 1-based positions """
        return ''.join('{}\t{}\t{}\n'.format(path, i + 1, d) for i, d in enumerate(depths))

    def awk(self, text):
        proc = subprocess.run(['awk', DEPTH_TO_BED_AWK], input=text,
                              capture_output=True, text=True, check=True)
        return proc.stdout.splitlines(keepends=True)

    def runs(self, depths, min_length=1, offsets=None, path='R#0#chr1'):
        return list(merge_sorted_bed_stream(self.awk(self.depth_lines(path, depths)),
                                            min_length, offsets))

    def test_all_zero_path(self):
        # a length-10 path with no other assembly aligned is exactly (0, 10), NOT (1, 11)
        self.assertEqual(self.runs([0] * 10), [('R#0#chr1', 0, 10)])

    def test_interior_gap(self):
        # 15bp reference, middle 5bp uncovered
        self.assertEqual(self.runs([1] * 5 + [0] * 5 + [1] * 5), [('R#0#chr1', 5, 10)])

    def test_run_touching_last_base(self):
        depths = [1] * 5 + [0] * 5
        runs = self.runs(depths)
        self.assertEqual(runs, [('R#0#chr1', 5, 10)])
        # end must equal the contig length, never length + 1
        self.assertEqual(runs[0][2], len(depths))

    def test_single_base_gap(self):
        self.assertEqual(self.runs([1, 0, 1]), [('R#0#chr1', 1, 2)])

    def test_min_length_filter(self):
        depths = [1] + [0] * 3 + [1] + [0] * 10
        self.assertEqual(self.runs(depths, min_length=5), [('R#0#chr1', 5, 15)])
        self.assertEqual(len(self.runs(depths, min_length=1)), 2)

    def test_multiple_paths(self):
        text = self.depth_lines('R#0#chr1', [0, 0]) + self.depth_lines('R#0#chr2', [0, 0, 0])
        self.assertEqual(list(merge_sorted_bed_stream(self.awk(text), 1)),
                         [('R#0#chr1', 0, 2), ('R#0#chr2', 0, 3)])

    def test_no_zero_depth(self):
        self.assertEqual(self.runs([1, 2, 3]), [])

    def test_vg_depth_positions_are_already_absolute(self):
        """ vg depth adds the subrange offset itself, so nothing downstream may add it again.
        this pins that the conversion is a pure 1-based -> 0-based shift """
        depths = ''.join('R#0#chr1\t{}\t0\n'.format(p) for p in range(1001, 1016))
        self.assertEqual(list(merge_sorted_bed_stream(self.awk(depths), 1)),
                         [('R#0#chr1', 1000, 1015)])

    def test_merge_is_streaming_over_adjacent_singletons(self):
        # merge_sorted_bed_stream must join abutting intervals, which is what awk emits per base
        lines = ['c\t0\t1\n', 'c\t1\t2\n', 'c\t2\t3\n', 'c\t10\t11\n']
        self.assertEqual(list(merge_sorted_bed_stream(lines, 1)), [('c', 0, 3), ('c', 10, 11)])


class TestSplitLog(unittest.TestCase):

    LOG = ('Assigned contig to chrXII: id=YPS128|chrXII  len=1027157 cov=0.972882 (vs 0.25) uf= infinity (vs 2)\n'
           'Assigned ref-contig to chrIII: id=S288C|chrIII  len=341580 cov=0.999701 (vs 0.5) uf= infinity (vs 2)\n'
           'Query contig is ambiguous: id=UWOPS034614|chrVII  len=632616 cov=0.407615 (vs 0.5) uf= infinity (vs 2)\n'
           ' Reference contig mappings:\n'
           '  chrVII: 353547\n')

    def test_parse(self):
        tmp_dir = tempfile.mkdtemp()
        try:
            path = os.path.join(tmp_dir, 'minigraph.split.log')
            with open(path, 'w') as log_file:
                log_file.write(self.LOG)
            assigned, ambiguous = parse_split_log(path)
            self.assertEqual(assigned, {('YPS128', 'chrXII', 0): 'chrXII',
                                        ('S288C', 'chrIII', 0): 'chrIII'})
            self.assertEqual(ambiguous, {('UWOPS034614', 'chrVII', 0)})
        finally:
            shutil.rmtree(tmp_dir)

    def test_fragment_names_resolve_to_contig_and_offset(self):
        """ the log names a fragment the way cactus named it internally; it has to land in the same
        (contig, offset) space as the graph and the input contig sizes, or a dropped fragment
        cannot be matched to its decision """
        tmp_dir = tempfile.mkdtemp()
        try:
            path = os.path.join(tmp_dir, 'minigraph.split.log')
            with open(path, 'w') as log_file:
                log_file.write('Assigned contig to chr6: id=HG002.1|chr6_sub_300000_597871  '
                               'len=297871 cov=0.93 (vs 0.5) uf= infinity (vs 2)\n')
                log_file.write('Query contig is ambiguous: id=HG002.1|chr1_sub_10_20_sub_2_5  '
                               'len=3 cov=0.1 (vs 0.5) uf=1.1 (vs 2)\n')
            assigned, ambiguous = parse_split_log(path)
            self.assertEqual(assigned, {('HG002.1', 'chr6', 300000): 'chr6'})
            # nested fragment names collapse to the summed offset, as hal2vg resolves them
            self.assertEqual(ambiguous, {('HG002.1', 'chr1', 12)})
        finally:
            shutil.rmtree(tmp_dir)

    def test_silently_dropped_contig_is_in_neither_bucket(self):
        """ a contig rgfa-split made no decision about appears in neither map.  conflating that
        silence with either bucket is exactly the bug this feature exists to surface """
        tmp_dir = tempfile.mkdtemp()
        try:
            path = os.path.join(tmp_dir, 'minigraph.split.log')
            with open(path, 'w') as log_file:
                log_file.write(self.LOG)
            assigned, ambiguous = parse_split_log(path)
            key = ('UWOPS034614', 'chrVI', 0)
            self.assertNotIn(key, assigned)
            self.assertNotIn(key, ambiguous)
        finally:
            shutil.rmtree(tmp_dir)


class TestComputeExclusions(unittest.TestCase):
    """ the classification, driven off small fixtures shaped like real coverage files """

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def write_cov(self, name, entries):
        path = os.path.join(self.tmp_dir, name)
        with open(path, 'w') as out_file:
            for pansn, contig, start, end in entries:
                out_file.write('{}\t{}\t{}\t{}\n'.format(pansn, contig, start, end))
        return path

    def read_bed(self, event, tag=''):
        """ tag selects the graph: '' is clip, '.full' the full graph, '.dN' the filtered one """
        path = os.path.join(self.tmp_dir, 'clipped' + tag + '.beds', '{}.bed'.format(event))
        with open(path, 'r') as in_file:
            return [tuple(l.rstrip('\n').split('\t')) for l in in_file]

    BASELINE = [('A', 'A#0', 'chr1', 0, 1000),
                ('A', 'A#0', 'chr2', 0, 500)]

    def test_clip_and_filter_are_attributed_separately(self):
        full = self.write_cov('c.full', [('A#0', 'chr1', 0, 1000)])
        clip = self.write_cov('c.clip', [('A#0', 'chr1', 100, 1000)])
        filt = self.write_cov('c.filter', [('A#0', 'chr1', 100, 900)])
        baseline = [('A', 'A#0', 'chr1', 0, 1000)]
        result = compute_exclusions(baseline, {'full': [full], 'clip': [clip], 'filter': [filt]},
                                    ['full', 'clip', 'filter'], None, ['chr1'], self.tmp_dir,
                                    filter_threshold=2)
        rows = self.read_bed('A', '.d2')
        self.assertIn(('chr1', '0', '100', 'clip'), rows)
        self.assertIn(('chr1', '900', '1000', 'filter'), rows)
        self.assertEqual(result['counts'][('A', 'clip')], [1, 100])
        self.assertEqual(result['counts'][('A', 'filter')], [1, 100])

    def test_absent_phase_produces_no_rows(self):
        """ the regression that matters most: under default options 'filter' never runs, and the
        report must not charge it the entire graph """
        full = self.write_cov('c.full', [('A#0', 'chr1', 0, 1000)])
        clip = self.write_cov('c.clip', [('A#0', 'chr1', 100, 1000)])
        baseline = [('A', 'A#0', 'chr1', 0, 1000)]
        result = compute_exclusions(baseline, {'full': [full], 'clip': [clip]}, ['full', 'clip'],
                                    None, ['chr1'], self.tmp_dir)
        self.assertNotIn(('A', 'filter'), result['counts'])
        self.assertEqual([r for r in self.read_bed('A') if r[3] == 'filter'], [])
        self.assertEqual(result['phases'], ['full', 'clip'])

    def test_contig_in_no_graph_is_dropped_not_clipped(self):
        full = self.write_cov('c.full', [('A#0', 'chr1', 0, 1000)])
        clip = self.write_cov('c.clip', [('A#0', 'chr1', 0, 1000)])
        result = compute_exclusions(self.BASELINE, {'full': [full], 'clip': [clip]},
                                    ['full', 'clip'], ({}, set()), ['chr1'], self.tmp_dir)
        rows = self.read_bed('A')
        self.assertIn(('chr2', '0', '500', 'unassigned'), rows)
        self.assertEqual(result['counts'][('A', 'unassigned')], [1, 500])

    def test_split_log_distinguishes_ambiguous_from_unassigned(self):
        full = self.write_cov('c.full', [('A#0', 'chr1', 0, 1000)])
        baseline = self.BASELINE + [('A', 'A#0', 'chr3', 0, 300)]
        split_log = ({}, {('A', 'chr2', 0)})
        result = compute_exclusions(baseline, {'full': [full]}, ['full'], split_log,
                                    ['chr1'], self.tmp_dir)
        rows = dict((r[0], r[3]) for r in self.read_bed('A', '.full'))
        self.assertEqual(rows['chr2'], 'ambiguous')
        self.assertEqual(rows['chr3'], 'unassigned')

    def test_assigned_to_missing_chromosome_graph(self):
        """ --refContigs without --otherContig drops whole chromosomes after assignment """
        full = self.write_cov('c.full', [('A#0', 'chr1', 0, 1000)])
        split_log = ({('A', 'chr2', 0): 'chr2'}, set())
        result = compute_exclusions(self.BASELINE, {'full': [full]}, ['full'], split_log,
                                    ['chr1'], self.tmp_dir)
        self.assertEqual(dict((r[0], r[3]) for r in self.read_bed('A', '.full'))['chr2'],
                         'no_chromosome_graph')

    def test_contig_folded_into_chrOther_is_unaligned_not_no_graph(self):
        """ the GRCh38 default: rgfa-split logs the fine-grained ref contig, cactus folds it into
        chrOther, which IS built.  a contig that reached chrOther and lost all its sequence there
        is unaligned, not no_chromosome_graph, and files under chrOther """
        full = self.write_cov('c.full', [('A#0', 'chr1', 0, 1000)])
        # chr2 was assigned to a random contig that has no graph of its own, folded into chrOther
        split_log = ({('A', 'chr2', 0): 'chr14_GL000009v2_random'}, set())
        result = compute_exclusions(self.BASELINE, {'full': [full]}, ['full'], split_log,
                                    ['chr1', 'chrOther'], self.tmp_dir, other_contig='chrOther')
        self.assertEqual(dict((r[0], r[3]) for r in self.read_bed('A', '.full'))['chr2'],
                         'unaligned')
        # and it files under the graph it reached, not the phantom random contig
        chrom_of = dict((contig, chrom) for chrom, _ev, contig, *_ in result['contig_rows'])
        self.assertEqual(chrom_of['chr2'], 'chrOther')
        # not exempt from the warning, unlike a genuine no_chromosome_graph
        self.assertIn(('A', 'unaligned'), result['counts'])

    def test_partially_covered_contig_is_never_split_labelled(self):
        """ any full-phase coverage proves the contig reached a graph, whatever the log says """
        full = self.write_cov('c.full', [('A#0', 'chr1', 0, 900)])
        baseline = [('A', 'A#0', 'chr1', 0, 1000)]
        split_log = ({}, {('A', 'chr1', 0)})
        compute_exclusions(baseline, {'full': [full]}, ['full'], split_log, ['chr1'], self.tmp_dir)
        self.assertEqual(dict((r[0], r[3]) for r in self.read_bed('A', '.full'))['chr1'], 'unaligned')

    def test_degraded_baseline_suppresses_dropped(self):
        full = self.write_cov('c.full', [('A#0', 'chr1', 0, 1000)])
        clip = self.write_cov('c.clip', [('A#0', 'chr1', 100, 1000)])
        result = compute_exclusions(None, {'full': [full], 'clip': [clip]}, ['full', 'clip'],
                                    None, ['chr1'], self.tmp_dir)
        self.assertFalse(result['report_dropped'])
        self.assertEqual(result['counts'][('A', 'clip')], [1, 100])
        self.assertNotIn(('A', 'unassigned'), result['counts'])

    def test_every_genome_gets_a_file_even_with_no_loss(self):
        full = self.write_cov('c.full', [('A#0', 'chr1', 0, 1000), ('A#0', 'chr2', 0, 500)])
        compute_exclusions(self.BASELINE, {'full': [full]}, ['full'], None, ['chr1'], self.tmp_dir)
        self.assertEqual(self.read_bed('A', '.full'), [])
        self.assertTrue(os.path.exists(os.path.join(self.tmp_dir, 'clipped.full.beds', 'A.bed')))

    def test_single_base_filter_fragments_are_enumerated(self):
        """ the allele-frequency filter removes private alleles, which are often 1bp.  they are
        real excluded sequence and must appear in the BED like anything else """
        full = self.write_cov('c.full', [('A#0', 'chr1', 0, 1000)])
        clip = self.write_cov('c.clip', [('A#0', 'chr1', 0, 1000)])
        filt = self.write_cov('c.filter', [('A#0', 'chr1', 0, 10), ('A#0', 'chr1', 11, 1000)])
        baseline = [('A', 'A#0', 'chr1', 0, 1000)]
        result = compute_exclusions(baseline, {'full': [full], 'clip': [clip], 'filter': [filt]},
                                    ['full', 'clip', 'filter'], None, ['chr1'], self.tmp_dir,
                                    filter_threshold=2)
        self.assertIn(('chr1', '10', '11', 'filter'), self.read_bed('A', '.d2'))
        self.assertEqual(result['counts'][('A', 'filter')], [1, 1])

    def test_each_graph_gets_its_own_cumulative_bed(self):
        """ the whole point of splitting by graph: what is missing from the clipped graph must not
        include what the later frequency filter took, and the filtered graph must include both """
        full = self.write_cov('c.full', [('A#0', 'chr1', 0, 1000)])
        clip = self.write_cov('c.clip', [('A#0', 'chr1', 100, 1000)])
        filt = self.write_cov('c.filter', [('A#0', 'chr1', 100, 900)])
        baseline = [('A', 'A#0', 'chr1', 0, 1000),
                    ('A', 'A#0', 'chr2', 0, 500)]
        compute_exclusions(baseline, {'full': [full], 'clip': [clip], 'filter': [filt]},
                           ['full', 'clip', 'filter'], ({}, set()), ['chr1'], self.tmp_dir,
                           filter_threshold=2)
        full_bed = self.read_bed('A', '.full')
        clip_bed = self.read_bed('A', '')
        filt_bed = self.read_bed('A', '.d2')
        # chr2 never reached a graph, so it is missing from all three
        self.assertIn(('chr2', '0', '500', 'unassigned'), full_bed)
        self.assertIn(('chr2', '0', '500', 'unassigned'), clip_bed)
        self.assertIn(('chr2', '0', '500', 'unassigned'), filt_bed)
        # the clip removal is absent from the full graph's set, present in the other two
        self.assertNotIn(('chr1', '0', '100', 'clip'), full_bed)
        self.assertIn(('chr1', '0', '100', 'clip'), clip_bed)
        self.assertIn(('chr1', '0', '100', 'clip'), filt_bed)
        # the filter removal appears only in the filtered graph's set
        self.assertNotIn(('chr1', '900', '1000', 'filter'), clip_bed)
        self.assertIn(('chr1', '900', '1000', 'filter'), filt_bed)
        self.assertTrue(len(full_bed) < len(clip_bed) < len(filt_bed))

    def test_coverage_outside_baseline_is_counted_not_crashed(self):
        full = self.write_cov('c.full', [('A#0', 'chr1', 0, 1200)])
        baseline = [('A', 'A#0', 'chr1', 0, 1000)]
        result = compute_exclusions(baseline, {'full': [full]}, ['full'], None,
                                    ['chr1'], self.tmp_dir)
        self.assertEqual(result['outside_baseline_bp'], 200)

    def test_graph_key_with_no_baseline_row_is_an_orphan_not_a_crash(self):
        full = self.write_cov('c.full', [('A#0', 'chr1', 0, 1000), ('Z#0', 'chrZ', 0, 77)])
        baseline = [('A', 'A#0', 'chr1', 0, 1000)]
        result = compute_exclusions(baseline, {'full': [full]}, ['full'], None,
                                    ['chr1'], self.tmp_dir)
        self.assertEqual(result['orphan_paths'], 1)
        self.assertEqual(result['orphan_bp'], 77)


class TestEventFilename(unittest.TestCase):

    def test_accepts_normal_names(self):
        self.assertEqual(safe_event_filename('HG002.1'), 'HG002.1')

    def test_rejects_path_separators(self):
        for bad in ['a/b', '..', '.']:
            with self.assertRaises(RuntimeError):
                safe_event_filename(bad)




class TestDetectProblems(unittest.TestCase):
    """ the WARNING heuristics.  the hard part is not spotting a bad run, it is staying quiet on a
    deliberately lossy one -- calibrated against real yeast and primate runs """

    @staticmethod
    def result(counts, input_bp, whole_dropped=None, outside=0, orphans=0):
        return {'counts': counts, 'input_bp': input_bp, 'whole_dropped': whole_dropped or {},
                'outside_baseline_bp': outside, 'orphan_paths': orphans, 'orphan_keys': []}

    def test_silent_on_a_healthy_panel(self):
        # every genome loses ~10% to the frequency filter, as the primates run does
        counts = dict(((s, 'filter'), [100, 1000000]) for s in ['a', 'b', 'c', 'd'])
        result = self.result(counts, dict((s, 10000000) for s in ['a', 'b', 'c', 'd']))
        self.assertEqual(detect_problems(result), [])

    def test_flags_a_genome_that_reached_no_graph(self):
        counts = {('bad', 'ambiguous'): [4, 3000000], ('ok', 'clip'): [2, 100000]}
        result = self.result(counts, {'bad': 10000000, 'ok': 10000000})
        problems = detect_problems(result)
        self.assertTrue(any('bad' in p and 'no graph at all' in p for p in problems))
        self.assertFalse(any(p.startswith('ok') for p in problems))

    def test_silent_on_uniform_deliberate_loss(self):
        """ --refContigs drops whole chromosomes from every genome on purpose.  that is
        no_chromosome_graph, it is uniform, and it must not raise anything """
        counts = dict(((s, 'no_chromosome_graph'), [3, 2400000]) for s in ['a', 'b', 'c', 'd'])
        result = self.result(counts, dict((s, 10000000) for s in ['a', 'b', 'c', 'd']))
        self.assertEqual(detect_problems(result), [])

    def test_flags_an_outlier_genome(self):
        counts = dict(((s, 'clip'), [10, 500000]) for s in ['a', 'b', 'c'])
        counts[('bad', 'clip')] = [10, 4000000]
        input_bp = dict((s, 10000000) for s in ['a', 'b', 'c', 'bad'])
        problems = detect_problems(self.result(counts, input_bp))
        self.assertTrue(any('bad' in p and 'differently from the rest' in p for p in problems))

    def test_names_whole_contigs_that_vanished(self):
        result = self.result({}, {'bad': 10000000},
                             whole_dropped={'bad': [('chrVII', 632616, 'ambiguous')]})
        problems = detect_problems(result)
        self.assertTrue(any('chrVII' in p and '632616' in p for p in problems))

    def test_flags_a_reference_nothing_aligns_to(self):
        result = self.result({}, {'ref': 10000000})
        problems = detect_problems(result, refgap_bp={'ref': (10, 5000000)},
                                   ref_input_bp=10000000)
        self.assertTrue(any('no other assembly aligned' in p for p in problems))

    def test_flags_internal_inconsistency(self):
        for kwargs in ({'outside': 1234}, {'orphans': 3}):
            problems = detect_problems(self.result({}, {'a': 1000}, **kwargs))
            self.assertTrue(any(p.startswith('INTERNAL') for p in problems))

    def test_no_division_by_zero_on_an_empty_genome(self):
        self.assertEqual(detect_problems(self.result({}, {'a': 0})), [])

if __name__ == '__main__':
    unittest.main()
