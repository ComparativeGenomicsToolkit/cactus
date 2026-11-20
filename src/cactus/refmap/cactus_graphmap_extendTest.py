#!/usr/bin/env python3

"""
Unit tests for the pure (non-Toil) logic behind extending a pangenome: the PanSN <-> cactus name
round trip on a published GAF, and splitting a merged GAF back into the per-genome pieces it was
concatenated from.

Both exist because minigraph GAF is in stable coordinates, which adding genomes to a graph does not
change, so a genome's existing mappings can be re-derived against the extended graph instead of
being recomputed.  That only holds if the published GAF round trips exactly, which is what these
pin down.

These are fast and offline.  The end-to-end extension is covered by evolverTest.py.

TestGafOverlapFilterConfig pins down which GAF overlap filter command each <graphmap> setting runs
(the exact filter, gaffilter -x, is on by default; "0" must give the stock command).

TestAddGenomeToSeqFile checks that the seqfile cactus-graphmap adds the minigraph assembly to keeps
the pangenome default branch length.

TestPruneOtherContigs, at the end, checks the prune of what a --mgSplit wholeGenomeRef construction
added to other chromosomes' reference (<graphmap_split pruneOtherContigs>): which segments and links
go, that the reference never does, that both namings and the --otherContig lump are handled, and
that anything but an rGFA is refused.
"""

import copy
import gzip
import os
import re
import shutil
import tempfile
import unittest
import xml.etree.ElementTree as ET

from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.progressive.progressive_decomposition import parse_seqfile
from cactus.refmap.cactus_graphmap import (
    gaf_from_pansn, gaf_to_pansn, pansn_to_event_map, split_gaf_file_by_event, trim_unstable_gaf,
    gaf_to_unstable_commands, overlap_filter_exact, check_overlap_filter_config, apply_mgsplit_filter_overrides,
    add_genome_to_seqfile)
from cactus.refmap.cactus_graphmap_split import (
    prune_other_contigs_gfa, prune_other_contigs_enabled, prune_other_contigs_memory, prune_other_contigs_plan,
    prune_other_contigs_report, prune_own_contigs, scan_rgfa, gfa_ref_contig, PRUNE_MIN_MEMORY)


# a stable GAF line as minigraph writes it: query name, then a path of stable segments with
# orientation marks, then the tags.  only columns 1 and 6 carry sequence names
def gaf_line(query, path, tags='60\ttp:A:P\tcm:i:100'):
    return '{}\t1000\t0\t900\t+\t{}\t2000\t10\t910\t880\t900\t{}\n'.format(query, path, tags)


class TestPansnToEventMap(unittest.TestCase):

    def test_haploid_and_diploid(self):
        self.assertEqual(pansn_to_event_map({'S288C', 'HG002#1', 'HG002#2'}),
                         {'S288C#0': 'S288C', 'HG002#1': 'HG002#1', 'HG002#2': 'HG002#2'})

    def test_explicit_hap_zero_beats_bare_name(self):
        # event_to_pansn_prefix maps both S288C and S288C#0 onto S288C#0.  whichever wins, the map
        # has to be single valued or the round trip would be ambiguous
        prefix_map = pansn_to_event_map({'S288C', 'S288C#0'})
        self.assertEqual(len(prefix_map), 1)
        self.assertIn(prefix_map['S288C#0'], ('S288C', 'S288C#0'))

    def test_dots_in_the_sample(self):
        # a '.' is never a haplotype in a genome name: the old SAMPLE.N spelling is only read in the seqfile
        self.assertEqual(pansn_to_event_map({'KOLF2.1J', 'GCA_054883195.1', 'HG002.1'}),
                         {'KOLF2.1J#0': 'KOLF2.1J', 'GCA_054883195.1#0': 'GCA_054883195.1', 'HG002.1#0': 'HG002.1'})

    def test_non_numeric_suffix_is_part_of_the_sample(self):
        # HG002.pat is not a haplotype suffix, so the whole thing is the sample
        self.assertEqual(pansn_to_event_map({'HG002.pat'}), {'HG002.pat#0': 'HG002.pat'})


class TestGafPansnRoundTrip(unittest.TestCase):

    def round_trip(self, names, cactus_gaf):
        """ cactus -> PanSN -> cactus, which is the path a reused GAF actually takes """
        with tempfile.TemporaryDirectory() as work_dir:
            in_path = os.path.join(work_dir, 'in.gaf')
            pansn_path = os.path.join(work_dir, 'pansn.gaf')
            back_path = os.path.join(work_dir, 'back.gaf')
            with open(in_path, 'w') as in_file:
                in_file.write(cactus_gaf)
            gaf_to_pansn(in_path, pansn_path)
            gaf_from_pansn(names, pansn_path, back_path)
            with open(pansn_path) as pansn_file, open(back_path) as back_file:
                return pansn_file.read(), back_file.read()

    def test_forward_and_back(self):
        names = {'S288C', 'SK1'}
        gaf = gaf_line('id=SK1|chrI', '>id=S288C|chrI>id=S288C|chrII')
        pansn, back = self.round_trip(names, gaf)
        self.assertIn('SK1#0#chrI', pansn)
        self.assertIn('>S288C#0#chrI>S288C#0#chrII', pansn)
        self.assertEqual(back, gaf)

    def test_reverse_steps(self):
        names = {'HG002#1', 'HG002#2', 'GRCh38'}
        gaf = gaf_line('id=HG002#2|chr1', '<id=GRCh38|chr1>id=HG002#1|chr1<id=GRCh38|chr1')
        pansn, back = self.round_trip(names, gaf)
        self.assertIn('<GRCh38#0#chr1>HG002#1#chr1<GRCh38#0#chr1', pansn)
        self.assertEqual(back, gaf)

    def test_stable_subrange_on_a_step(self):
        # a step can name a sub-interval of the stable sequence; only the name part is rewritten
        names = {'S288C', 'SK1'}
        gaf = gaf_line('id=SK1|chrI', '>id=S288C|chrI:100-2000')
        pansn, back = self.round_trip(names, gaf)
        self.assertIn('>S288C#0#chrI:100-2000', pansn)
        self.assertEqual(back, gaf)

    def test_contig_name_containing_a_hash(self):
        # PanSN phase blocks put a fourth '#' field on a path name.  going back, only the first two
        # '#' belong to the prefix
        names = {'HG002#1', 'GRCh38'}
        gaf = gaf_line('id=HG002#1|chr1#0', '>id=GRCh38|chr1#0')
        pansn, back = self.round_trip(names, gaf)
        self.assertIn('HG002#1#chr1#0', pansn)
        self.assertEqual(back, gaf)

    def test_dots_in_the_sample(self):
        names = {'KOLF2.1J', 'GCA_054883195.1'}
        gaf = gaf_line('id=KOLF2.1J|chr1', '>id=GCA_054883195.1|chr1<id=KOLF2.1J|chr1')
        pansn, back = self.round_trip(names, gaf)
        self.assertIn('KOLF2.1J#0#chr1', pansn)
        self.assertIn('>GCA_054883195.1#0#chr1<KOLF2.1J#0#chr1', pansn)
        self.assertEqual(back, gaf)

    def from_pansn(self, names, gaf):
        with tempfile.TemporaryDirectory() as work_dir:
            in_path = os.path.join(work_dir, 'in.gaf')
            out_path = os.path.join(work_dir, 'out.gaf')
            with open(in_path, 'w') as in_file:
                in_file.write(gaf)
            gaf_from_pansn(names, in_path, out_path)
            with open(out_path) as out_file:
                return out_file.read()

    def test_genome_spelled_with_its_hap_zero_is_updated(self):
        # a cactus-named GAF from a cactus that kept the #0 of a genome given as S288C#0, read against a
        # seqfile whose names are normalized to S288C
        gaf = gaf_line('id=SK1|chrI', '>id=S288C#0|chrI<id=SK1|chrI')
        self.assertEqual(self.from_pansn({'S288C', 'SK1'}, gaf), gaf_line('id=SK1|chrI', '>id=S288C|chrI<id=SK1|chrI'))

    def test_dot_is_not_a_haplotype(self):
        # the old SAMPLE.N spelling is only read in a seqfile, so in a GAF it is a genome of that name
        gaf = gaf_line('id=HG002.1|chr1', '>id=GRCh38|chr1<id=HG002.1|chr1')
        self.assertEqual(self.from_pansn({'HG002#1', 'GRCh38'}, gaf), gaf)

    def test_contig_name_containing_a_pipe(self):
        names = {'S288C', 'SK1'}
        gaf = gaf_line('id=SK1|ctg|1', '>id=S288C|chrI')
        pansn, back = self.round_trip(names, gaf)
        self.assertEqual(back, gaf)

    def test_already_cactus_named_gaf_passes_through(self):
        # an older cactus published the GAF in its own naming; reading one back must not mangle it
        names = {'S288C', 'SK1'}
        gaf = gaf_line('id=SK1|chrI', '>id=S288C|chrI')
        with tempfile.TemporaryDirectory() as work_dir:
            in_path = os.path.join(work_dir, 'in.gaf')
            out_path = os.path.join(work_dir, 'out.gaf')
            with open(in_path, 'w') as in_file:
                in_file.write(gaf)
            gaf_from_pansn(names, in_path, out_path)
            with open(out_path) as out_file:
                self.assertEqual(out_file.read(), gaf)

    def test_unknown_prefix_is_left_alone(self):
        # a contig that merely looks like a PanSN prefix belongs to no genome and must not be
        # rewritten into a name that exists in no graph
        gaf = gaf_line('weird#name#ctg', '>other#thing#ctg')
        with tempfile.TemporaryDirectory() as work_dir:
            in_path = os.path.join(work_dir, 'in.gaf')
            out_path = os.path.join(work_dir, 'out.gaf')
            with open(in_path, 'w') as in_file:
                in_file.write(gaf)
            gaf_from_pansn({'S288C'}, in_path, out_path)
            with open(out_path) as out_file:
                self.assertEqual(out_file.read(), gaf)


class TestSplitGafByEvent(unittest.TestCase):

    def split(self, gaf_text, names, gzipped=False):
        work_dir = tempfile.mkdtemp()
        gaf_path = os.path.join(work_dir, 'merged.gaf.gz' if gzipped else 'merged.gaf')
        if gzipped:
            with gzip.open(gaf_path, 'wt') as gaf_file:
                gaf_file.write(gaf_text)
        else:
            with open(gaf_path, 'w') as gaf_file:
                gaf_file.write(gaf_text)
        shard_dir = os.path.join(work_dir, 'shards')
        os.makedirs(shard_dir)
        shard_paths, dropped = split_gaf_file_by_event(gaf_path, names, shard_dir)
        contents = {}
        for event, shard_path in shard_paths.items():
            with open(shard_path) as shard_file:
                contents[event] = shard_file.read()
        return contents, dropped

    def test_partition_is_exact(self):
        # every line lands in exactly one shard, in order: this is what makes the re-derived PAF
        # identical to the original when nothing has been added
        names = {'S288C', 'SK1', 'Y12'}
        lines = [gaf_line('S288C#0#chrI', '>S288C#0#chrI'),
                 gaf_line('S288C#0#chrII', '>S288C#0#chrII'),
                 gaf_line('SK1#0#chrI', '>S288C#0#chrI'),
                 gaf_line('Y12#0#chrI', '>S288C#0#chrI')]
        contents, dropped = self.split(''.join(lines), names)
        self.assertEqual(set(contents), names)
        self.assertEqual(dropped, set())
        self.assertEqual(contents['S288C'], lines[0] + lines[1])
        self.assertEqual(contents['SK1'], lines[2])
        self.assertEqual(''.join(contents[e] for e in ['S288C', 'SK1', 'Y12']), ''.join(lines))

    def test_gzipped_input(self):
        names = {'S288C'}
        line = gaf_line('S288C#0#chrI', '>S288C#0#chrI')
        contents, dropped = self.split(line, names, gzipped=True)
        self.assertEqual(contents, {'S288C': line})

    def test_interleaved_records_still_group(self):
        # append mode means the grouping does not depend on the concatenation order
        names = {'S288C', 'SK1'}
        a1 = gaf_line('S288C#0#chrI', '>S288C#0#chrI')
        b1 = gaf_line('SK1#0#chrI', '>S288C#0#chrI')
        a2 = gaf_line('S288C#0#chrII', '>S288C#0#chrII')
        contents, _ = self.split(a1 + b1 + a2, names)
        self.assertEqual(contents['S288C'], a1 + a2)
        self.assertEqual(contents['SK1'], b1)

    def test_diploid_haplotypes_are_separate_shards(self):
        names = {'HG002#1', 'HG002#2'}
        h1 = gaf_line('HG002#1#chr1', '>GRCh38#0#chr1')
        h2 = gaf_line('HG002#2#chr1', '>GRCh38#0#chr1')
        contents, _ = self.split(h1 + h2, names)
        self.assertEqual(contents, {'HG002#1': h1, 'HG002#2': h2})

    def test_hap_zero_spelling_finds_its_genome(self):
        names = {'S288C', 'SK1'}
        a = gaf_line('id=S288C#0|chrI', '>id=S288C#0|chrI')
        b = gaf_line('id=SK1|chrI', '>id=S288C#0|chrI')
        contents, dropped = self.split(a + b, names)
        self.assertEqual(contents, {'S288C': a, 'SK1': b})
        self.assertEqual(dropped, set())

    def test_dot_is_not_a_haplotype(self):
        # id=HG002.1| is a genome named HG002.1, which is not in the seqfile
        names = {'HG002#1', 'HG002#2'}
        h1 = gaf_line('id=HG002.1|chr1', '>id=GRCh38|chr1')
        h2 = gaf_line('id=HG002#2|chr1', '>id=GRCh38|chr1')
        contents, dropped = self.split(h1 + h2, names)
        self.assertEqual(contents, {'HG002#2': h2})
        self.assertEqual(dropped, {'HG002.1'})

    def test_genome_not_in_seqfile_is_dropped_and_reported(self):
        names = {'S288C'}
        keep = gaf_line('S288C#0#chrI', '>S288C#0#chrI')
        drop = gaf_line('GONE#0#chrI', '>S288C#0#chrI')
        contents, dropped = self.split(keep + drop, names)
        self.assertEqual(contents, {'S288C': keep})
        self.assertEqual(dropped, {'GONE'})

    def test_cactus_named_gaf_splits_too(self):
        names = {'S288C', 'SK1'}
        a = gaf_line('id=S288C|chrI', '>id=S288C|chrI')
        b = gaf_line('id=SK1|chrI', '>id=S288C|chrI')
        contents, dropped = self.split(a + b, names)
        self.assertEqual(contents, {'S288C': a, 'SK1': b})
        self.assertEqual(dropped, set())


class TestTrimUnstableGaf(unittest.TestCase):
    """ gaf2paf reads a record's path start as an offset into its first step.  Extending a graph
    splits nodes, so a reused mapping's offset can end up past the first of the finer nodes that
    replaced the one it was made against -- which gaf2paf asserts on.  Trimming the steps that hold
    none of the alignment puts the offsets back where gaf2paf expects them. """

    NODE_LENS = {'s1': 100, 's2': 50, 's3': 200, 's4': 30, 's5': 80}

    def trim(self, records):
        """ run trim_unstable_gaf over (path, path_len, path_start, path_end) tuples, returning the
        same tuples back """
        with tempfile.TemporaryDirectory() as work_dir:
            lengths_path = os.path.join(work_dir, 'lens.tsv')
            with open(lengths_path, 'w') as lengths_file:
                for node, length in self.NODE_LENS.items():
                    lengths_file.write('{}\t{}\n'.format(node, length))
            in_path = os.path.join(work_dir, 'in.gaf')
            out_path = os.path.join(work_dir, 'out.gaf')
            with open(in_path, 'w') as in_file:
                for path, path_len, path_start, path_end in records:
                    in_file.write('q\t1000\t0\t900\t+\t{}\t{}\t{}\t{}\t880\t900\t60\tcg:Z:900M\n'.format(
                        path, path_len, path_start, path_end))
            trimmed = trim_unstable_gaf(in_path, out_path, lengths_path)
            out = []
            with open(out_path) as out_file:
                for line in out_file:
                    toks = line.rstrip('\n').split('\t')
                    out.append((toks[5], int(toks[6]), int(toks[7]), int(toks[8])))
            return out, trimmed

    def assert_gaf2paf_invariant(self, record):
        """ what gaf2paf assumes: the start offset is inside the first step and the end offset is
        inside the last, and the steps still add up to the stated path length """
        path, path_len, path_start, path_end = record
        steps = [s[1:] for s in re.findall(r'[<>][^<>]+', path)]
        node_lens = [self.NODE_LENS[s] for s in steps]
        self.assertEqual(sum(node_lens), path_len)
        self.assertLess(path_start, node_lens[0])
        self.assertGreater(path_end, path_len - node_lens[-1])

    def test_nothing_to_trim_passes_through(self):
        # the mapping was made against this very graph: every offset is already in its end step
        record = ('>s1>s2>s3', 350, 40, 300)
        out, trimmed = self.trim([record])
        self.assertEqual(out, [record])
        self.assertEqual(trimmed, 0)

    def test_offset_past_the_first_step(self):
        # s1+s2 replaced one 150bp node, so an offset of 120 into it now lands in s2
        out, trimmed = self.trim([('>s1>s2>s3', 350, 120, 300)])
        self.assertEqual(trimmed, 1)
        self.assertEqual(out, [('>s2>s3', 250, 20, 200)])
        self.assert_gaf2paf_invariant(out[0])

    def test_alignment_ends_before_the_last_steps(self):
        out, trimmed = self.trim([('>s3>s2>s1', 350, 10, 200)])
        self.assertEqual(trimmed, 1)
        self.assertEqual(out, [('>s3', 200, 10, 200)])
        self.assert_gaf2paf_invariant(out[0])

    def test_trims_both_ends(self):
        out, trimmed = self.trim([('>s1>s2>s3>s4>s5', 460, 150, 350)])
        self.assertEqual(trimmed, 1)
        # s1 and s2 are wholly before the start, s4 and s5 wholly after the end
        self.assertEqual(out, [('>s3', 200, 0, 200)])
        self.assert_gaf2paf_invariant(out[0])

    def test_offsets_keep_their_span(self):
        # trimming moves the window, it must never resize it
        for record in [('>s1>s2>s3', 350, 120, 300), ('>s1>s2>s3>s4>s5', 460, 150, 350)]:
            out, _ = self.trim([record])
            self.assertEqual(out[0][3] - out[0][2], record[3] - record[2])

    def test_reverse_steps_trim_the_same_way(self):
        # the offsets run along the path as written, so orientation does not enter into it
        out, trimmed = self.trim([('<s1<s2<s3', 350, 120, 300)])
        self.assertEqual(trimmed, 1)
        self.assertEqual(out, [('<s2<s3', 250, 20, 200)])

    def test_single_step_is_never_emptied(self):
        out, trimmed = self.trim([('>s1', 100, 100, 100)])
        self.assertEqual(out, [('>s1', 100, 100, 100)])
        self.assertEqual(trimmed, 0)

    def test_non_path_lines_pass_through(self):
        # a stable path (no orientation marks) is not ours to touch
        with tempfile.TemporaryDirectory() as work_dir:
            lengths_path = os.path.join(work_dir, 'lens.tsv')
            open(lengths_path, 'w').write('s1\t100\n')
            in_path = os.path.join(work_dir, 'in.gaf')
            out_path = os.path.join(work_dir, 'out.gaf')
            text = gaf_line('q', 'chr1')
            open(in_path, 'w').write(text)
            self.assertEqual(trim_unstable_gaf(in_path, out_path, lengths_path), 0)
            with open(out_path) as out_file:
                self.assertEqual(out_file.read(), text)


class TestGafOverlapFilterConfig(unittest.TestCase):
    """ the gaf2unstable | gaffilter command stable_gaf_to_paf() runs, per <graphmap> setting """

    # what the default <graphmap> (minMAPQ 5, minGAFBlockLength 250000, minIdentity 0.5, ratio 5,
    # min length ratio 0.25) gives with the exact filter off: the command every release before it ran
    STOCK = [['gaf2unstable', 'g.gaf', '-g', 'g.gfa', '-o', 'g.len'],
             ['gaffilter', '-', '-r', '5.0', '-m', '0.25', '-q', '5', '-b', '250000', '-i', '0.5']]

    def config(self, **attribs):
        root = ET.parse(os.path.join(cactusRootPath(), 'cactus_progressive_config.xml')).getroot()
        root.find('graphmap').attrib.update(attribs)
        return root

    def commands(self, root, reference='S288C'):
        return gaf_to_unstable_commands(root, 'g.gaf', 'g.gfa', 'g.len', reference)

    def test_exact_is_the_default(self):
        root = self.config()
        self.assertTrue(overlap_filter_exact(root))
        self.assertTrue(check_overlap_filter_config(root))
        cmd, junctions_path, summary_path = self.commands(root)
        self.assertEqual(cmd[0], self.STOCK[0] + ['-n', 'g.gfa.node_table.tsv'])
        self.assertEqual(cmd[1], self.STOCK[1] + ['-x', '--exact-nodes', 'g.gfa.node_table.tsv', '--exact-ref', 'S288C',
                                                  '--a0-paf-ratio', '5.0', '--a0-paf-min-overlap', '0.0', '--gap', '31000',
                                                  '--win-ident', '0.01', '--junctions', 'g.gaf.junctions.tsv',
                                                  '--exact-summary', 'g.gaf.exact.log'])
        self.assertEqual((junctions_path, summary_path), ('g.gaf.junctions.tsv', 'g.gaf.exact.log'))

    def test_off_is_stock(self):
        # "0", and a config from before the attribute existed, both run the stock filter, unchanged
        root = self.config(GAFOverlapFilterExact='0')
        self.assertFalse(check_overlap_filter_config(root))
        self.assertEqual(self.commands(root), (self.STOCK, None, None))
        del root.find('graphmap').attrib['GAFOverlapFilterExact']
        self.assertEqual(self.commands(root), (self.STOCK, None, None))

    def test_no_overlap_filter(self):
        # with GAFOverlapFilterRatio 0 there is no overlap filter to make exact: no error, no gaffilter
        root = self.config(GAFOverlapFilterRatio='0')
        self.assertFalse(overlap_filter_exact(root))
        self.assertEqual(self.commands(root), (self.STOCK[0], None, None))
        root = self.config()
        apply_mgsplit_filter_overrides(root)
        self.assertFalse(overlap_filter_exact(root))
        self.assertEqual(self.commands(root)[1:], (None, None))

    def test_exact_options(self):
        root = self.config(GAFOverlapFilterExactRatio='2', GAFOverlapFilterExactGuard='1', GAFOverlapFilterExactGap='21000',
                           GAFOverlapFilterExactOptions='--min-remainder 5000', collapse='all')
        gaffilter = self.commands(root, reference=None)[0][1]
        self.assertNotIn('--exact-ref', gaffilter)
        # collapsing skips filter_paf's -p stage, which the exact filter's notion of "new" models
        self.assertEqual(gaffilter[gaffilter.index('--a0-paf-ratio') + 1], '0')
        for opts in (['--exact-ratio', '2'], ['--guard'], ['--gap', '21000'], ['--min-remainder', '5000']):
            i = gaffilter.index(opts[0])
            self.assertEqual(gaffilter[i:i + len(opts)], opts)

    def test_win_ident(self):
        # on at 0.01 by default; any other margin is passed on as written
        gaffilter = self.commands(self.config(GAFOverlapFilterExactWinIdent='0.05'))[0][1]
        i = gaffilter.index('--win-ident')
        self.assertEqual(gaffilter[i:i + 2], ['--win-ident', '0.05'])
        # 0, empty, and a config from before the attribute existed all leave it off
        for value in ('0', '', None):
            root = self.config(GAFOverlapFilterExactWinIdent=value or '')
            if value is None:
                del root.find('graphmap').attrib['GAFOverlapFilterExactWinIdent']
            self.assertTrue(check_overlap_filter_config(root))
            self.assertNotIn('--win-ident', self.commands(root)[0][1])
        for value in ('-0.01', 'x'):
            with self.assertRaisesRegex(RuntimeError, 'GAFOverlapFilterExactWinIdent> must be a number'):
                check_overlap_filter_config(self.config(GAFOverlapFilterExactWinIdent=value))
        # it is a gaffilter -x option, so the stock filter neither checks nor passes it
        root = self.config(GAFOverlapFilterExact='0', GAFOverlapFilterExactWinIdent='x')
        self.assertFalse(check_overlap_filter_config(root))
        self.assertEqual(self.commands(root), (self.STOCK, None, None))
        with self.assertRaises(RuntimeError):
            check_overlap_filter_config(self.config(GAFOverlapFilterExactRatio='0'))

    def test_trim_removed(self):
        # gaffilter -t is gone: a config that still asks for it fails rather than silently getting -x
        for value in ('1', 'true'):
            with self.assertRaisesRegex(RuntimeError, 'GAFOverlapFilterTrim> has been removed'):
                check_overlap_filter_config(self.config(GAFOverlapFilterTrim=value))
            with self.assertRaises(RuntimeError):
                self.commands(self.config(GAFOverlapFilterTrim=value, GAFOverlapFilterExact='0'))
        # turned off (or its tuning attributes left behind) is harmless
        root = self.config(GAFOverlapFilterTrim='0', GAFOverlapFilterTrimEdge='5000', GAFOverlapFilterTrimMinMAPQ='20',
                           GAFOverlapFilterExact='0')
        self.assertEqual(self.commands(root), (self.STOCK, None, None))


class TestAddGenomeToSeqFile(unittest.TestCase):
    """ cactus-graphmap adds the minigraph assembly to the seqfile cactus-align then reads (the
    chrom-graphmap/*.gm.seqfile of --mgSplit).  the tree cactus-align gets must have the pangenome
    default branch length, which picks <caf annealingRounds>, rather than SeqFile's own 1.0 """

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        config_root = ET.parse(os.path.join(cactusRootPath(), 'cactus_progressive_config.xml')).getroot()
        # the seqfiles below use the old HG002.N spelling of a haplotype
        config_root.find('graphmap').set('dotHaplotypes', '1')
        self.config = ConfigWrapper(config_root)
        self.pg_len = self.config.getDefaultBranchLen(pangenome=True)

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def write_seqfile(self, lines):
        path = os.path.join(self.temp_dir, 'seqfile.txt')
        with open(path, 'w') as seq_file:
            seq_file.write('\n'.join(lines) + '\n')
        return path

    def branch_lengths(self, path):
        mc_tree = parse_seqfile(path, self.config, pangenome=True)[0]
        return {mc_tree.getName(leaf) : mc_tree.getWeight(mc_tree.getParent(leaf), leaf) for leaf in mc_tree.getLeaves()}

    def test_no_tree(self):
        # written back without a tree, so cactus-align makes its own star tree
        self.assertNotEqual(self.pg_len, 1.0)
        path = self.write_seqfile(['CHM13\tchm13.fa', 'HG002.1\tpat.fa', '*HG002#2\tmat.fa'])
        add_genome_to_seqfile(path, 'mg.fa', '_MINIGRAPH_', self.pg_len)
        with open(path) as seq_file:
            self.assertEqual(seq_file.read(), 'CHM13\tchm13.fa\nHG002.1\tpat.fa\n*HG002#2\tmat.fa\n_MINIGRAPH_\tmg.fa\n')
        self.assertEqual(self.branch_lengths(path), {name : self.pg_len for name in ['CHM13', 'HG002#1', 'HG002#2', '_MINIGRAPH_']})

    def test_tree(self):
        # a given tree keeps its lengths, and the minigraph assembly and a missing branch get the pangenome default
        path = self.write_seqfile(['(CHM13:0.01,HG002.1,HG002.2:0.03);', 'CHM13\tchm13.fa', 'HG002.1\tpat.fa', 'HG002.2\tmat.fa'])
        add_genome_to_seqfile(path, 'mg.fa', '_MINIGRAPH_', self.pg_len)
        self.assertEqual(self.branch_lengths(path), {'CHM13' : 0.01, 'HG002#1' : self.pg_len, 'HG002#2' : 0.03, '_MINIGRAPH_' : self.pg_len})
        # adding it again changes nothing
        with open(path) as seq_file:
            before = seq_file.read()
        add_genome_to_seqfile(path, 'mg.fa', '_MINIGRAPH_', self.pg_len)
        with open(path) as seq_file:
            self.assertEqual(seq_file.read(), before)



def rgfa_segment(name, seq, sn, sr, so=0):
    """ an rGFA S line as minigraph writes it """
    return 'S\t{}\t{}\tLN:i:{}\tSN:Z:{}\tSO:i:{}\tSR:i:{}\n'.format(name, seq, len(seq), sn, so, sr)

def rgfa_link(a, b, sr):
    return 'L\t{}\t+\t{}\t+\t0M\tSR:i:{}\n'.format(a, b, sr)


class TestPruneOtherContigs(unittest.TestCase):
    """ the prune of a --mgSplit wholeGenomeRef graph, built for one chromosome but against the whole
    reference, back to what that chromosome's construction could have added to its own reference
    (prune_other_contigs_plan) """

    # CHM13's chr4, chr10 and chr22, each as linear rank-0 segments, plus what chr4's haplotypes added:
    #   s9          a bubble on chr4                                  kept
    #   s10-s11     a component on chr10 alone (eg D4Z4 units)        dropped, with its 3 links
    #   s12         bridging chr4's s3 and chr22's s7                 kept, its link onto s7 dropped
    #   s13         attached to nothing                               dropped
    #   s14         a component on chr22 alone                        dropped, with its link
    # and three sample links between rank-0 segments: a deletion on chr4 (kept), a deletion on chr10
    # and a translocation from chr4 onto chr10 (both dropped)
    SEGMENTS = [('s1', 'ACGTACGTAC', 'chr4', 0), ('s2', 'ACGTACGTAC', 'chr4', 0), ('s3', 'ACGTACGTAC', 'chr4', 0),
                ('s4', 'GGGGCCCCAA', 'chr10', 0), ('s5', 'GGGGCCCCAA', 'chr10', 0), ('s6', 'GGGGCCCCAA', 'chr10', 0),
                ('s7', 'TTTTAAAA', 'chr22', 0), ('s8', 'TTTTAAAA', 'chr22', 0),
                ('s9', 'AAA', 'HG1#1|ctg1', 1),
                ('s10', 'CCCCCCCCCCCC', 'HG1#1|ctg1', 1), ('s11', 'CCCCCCCC', 'HG1#1|ctg1', 1),
                ('s12', 'GGGGG', 'HG2#1|ctg7', 2),
                ('s13', 'TT', 'HG2#1|ctg9', 2),
                ('s14', 'AAAAAAA', 'HG1#1|ctg3', 1)]
    LINKS = [('s1', 's2', 0), ('s2', 's3', 0), ('s4', 's5', 0), ('s5', 's6', 0), ('s7', 's8', 0),
             ('s1', 's9', 1), ('s9', 's2', 1),
             ('s4', 's10', 1), ('s10', 's11', 1), ('s11', 's5', 1),
             ('s3', 's12', 2), ('s12', 's7', 2),
             ('s8', 's14', 1),
             ('s1', 's3', 1),
             ('s4', 's6', 2),
             ('s3', 's4', 2)]
    KEPT_SEGMENTS = {'s1', 's2', 's3', 's4', 's5', 's6', 's7', 's8', 's9', 's12'}
    KEPT_LINKS = {('s1', 's2'), ('s2', 's3'), ('s4', 's5'), ('s5', 's6'), ('s7', 's8'), ('s1', 's9'), ('s9', 's2'),
                  ('s3', 's12'), ('s1', 's3')}

    def setUp(self):
        self.work_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.work_dir)

    def write_gfa(self, segments, links, name='g.gfa', pansn=False, gzipped=False):
        """ the toy graph, named as the cactus pipeline names it (id=CHM13|chr4) or in PanSN (CHM13#0#chr4) """
        def sn(name, rank):
            if rank == 0:
                return 'CHM13#0#' + name if pansn else 'id=CHM13|' + name
            return name.replace('|', '#') if pansn else 'id=' + name
        text = 'H\tVN:Z:1.0\n'
        text += ''.join(rgfa_segment(n, seq, sn(c, r), r) for n, seq, c, r in segments)
        text += ''.join(rgfa_link(a, b, r) for a, b, r in links)
        path = os.path.join(self.work_dir, name + ('.gz' if gzipped else ''))
        with (gzip.open(path, 'wt') if gzipped else open(path, 'w')) as gfa_file:
            gfa_file.write(text)
        return path

    def read_gfa(self, path):
        """ (S names, L (from, to) pairs, every line) of a GFA """
        with open(path) as gfa_file:
            lines = gfa_file.readlines()
        segs = set(line.split('\t')[1] for line in lines if line.startswith('S\t'))
        links = set(tuple(line.split('\t')[1:4:2]) for line in lines if line.startswith('L\t'))
        return segs, links, lines

    def prune(self, chrom, segments=None, links=None, bin_ref_contigs=None, **kwargs):
        in_path = self.write_gfa(segments or self.SEGMENTS, links or self.LINKS, **kwargs)
        out_path = os.path.join(self.work_dir, 'pruned.gfa')
        summary = prune_other_contigs_gfa(in_path, out_path, chrom, bin_ref_contigs)
        return summary, self.read_gfa(out_path)

    def assert_reference_kept(self, segments, links, kept_segs, kept_links):
        # every rank-0 segment and rank-0 link survives, whatever else goes
        self.assertTrue(set(n for n, _, _, r in segments if r == 0) <= kept_segs)
        self.assertTrue(set((a, b) for a, b, r in links if r == 0) <= kept_links)

    def test_other_chromosome_material_goes(self):
        summary, (segs, links, lines) = self.prune('chr4')
        self.assertEqual(segs, self.KEPT_SEGMENTS)
        self.assertEqual(links, self.KEPT_LINKS)
        self.assert_reference_kept(self.SEGMENTS, self.LINKS, segs, links)
        # anything that is neither S nor L is copied, and what is kept is copied unchanged and in order
        self.assertEqual(lines[0], 'H\tVN:Z:1.0\n')
        with open(os.path.join(self.work_dir, 'g.gfa')) as gfa_file:
            self.assertEqual([line for line in gfa_file if line in lines], lines)
        self.assertEqual(summary['own_contigs'], ['chr4'])
        self.assertEqual(summary['dropped_components'], {'chr10': (1, 2, 20), 'chr22': (1, 1, 7), 'none': (1, 1, 2)})
        self.assertEqual(summary['bridging_components'], {'chr22,chr4': (1, 1, 5)})
        self.assertEqual(summary['links_to_dropped_segments'], 4)
        self.assertEqual(summary['sample_links_onto_other_contigs'], {'chr10': 2, 'chr22': 1})
        self.assertEqual((summary['dropped_segments'], summary['dropped_bp'], summary['dropped_links']), (4, 29, 7))
        report = prune_other_contigs_report(summary)
        self.assertIn('dropped 1 non-reference component(s) attached only to chr10: 2 segments, 20 bp', report)
        self.assertIn('kept 1 bridging component(s) attached to chr22,chr4: 1 segments, 5 bp (their sample links onto '
                      'other contigs are dropped)', report)
        self.assertIn('dropped 1 non-reference component(s) attached to no reference contig: 1 segments, 2 bp', report)
        self.assertEqual(report[-1], 'kept 10 of 14 segments (29 bp dropped) and 9 of 16 links')

    def test_from_the_other_side(self):
        # chr10's graph keeps the chr10 component and the chr10 deletion, and loses chr4's bubble, the
        # bridge (which touches no chr10 reference), the chr4 deletion and the translocation, whose
        # other end is on chr4's reference
        summary, (segs, links, _) = self.prune('chr10')
        self.assertEqual(segs, {'s1', 's2', 's3', 's4', 's5', 's6', 's7', 's8', 's10', 's11'})
        self.assertEqual(links, {('s1', 's2'), ('s2', 's3'), ('s4', 's5'), ('s5', 's6'), ('s7', 's8'),
                                 ('s4', 's10'), ('s10', 's11'), ('s11', 's5'), ('s4', 's6')})
        self.assertEqual(summary['bridging_components'], {})
        self.assertEqual(summary['sample_links_onto_other_contigs'], {'chr4': 2})

    def test_pansn_and_gzip_prune_the_same(self):
        _, cactus = self.prune('chr4')
        for kwargs in ({'pansn': True}, {'gzipped': True}, {'pansn': True, 'gzipped': True}):
            summary, pruned = self.prune('chr4', **kwargs)
            self.assertEqual(pruned[:2], cactus[:2])
            self.assertEqual(summary['dropped_links'], 7)

    def test_contig_names(self):
        self.assertEqual(gfa_ref_contig('id=CHM13|chr4'), 'chr4')
        self.assertEqual(gfa_ref_contig('CHM13#0#chr4'), 'chr4')
        self.assertEqual(gfa_ref_contig('S288C#chrIV'), 'chrIV')
        self.assertEqual(gfa_ref_contig('chr4'), 'chr4')
        self.assertEqual(gfa_ref_contig('id=CHM13|chrUn#1'), 'chrUn#1')

    def test_other_contig_lump(self):
        # the --otherContig bin (chrOther) owns the contigs its own reference fasta holds: chrUn_a and
        # chrUn_b, not chr1.  A sample link between its two contigs is its own, and stays
        segments = [('a1', 'ACGT', 'chrUn_a', 0), ('a2', 'ACGT', 'chrUn_a', 0), ('b1', 'GGCC', 'chrUn_b', 0),
                    ('c1', 'TTTTTTTT', 'chr1', 0), ('c2', 'TTTTTTTT', 'chr1', 0),
                    ('x1', 'AA', 'HG1#1|u1', 1), ('x2', 'CCC', 'HG1#1|u2', 1)]
        links = [('a1', 'a2', 0), ('c1', 'c2', 0), ('a1', 'x1', 1), ('x1', 'a2', 1), ('c1', 'x2', 1), ('x2', 'c2', 1),
                 ('a2', 'b1', 1), ('b1', 'c1', 1)]
        calls = []
        def bin_ref_contigs():
            calls.append(1)
            return {'id=CHM13|chrUn_a', 'id=CHM13|chrUn_b', 'id=CHM13|chrUn_gone'}
        summary, (segs, kept_links, _) = self.prune('chrOther', segments, links, bin_ref_contigs)
        self.assertEqual(calls, [1])
        self.assertEqual(summary['own_contigs'], ['chrUn_a', 'chrUn_b'])
        self.assertEqual(segs, {'a1', 'a2', 'b1', 'c1', 'c2', 'x1'})
        self.assertEqual(kept_links, {('a1', 'a2'), ('c1', 'c2'), ('a1', 'x1'), ('x1', 'a2'), ('a2', 'b1')})
        self.assert_reference_kept(segments, links, segs, kept_links)
        # a plain set does as well as the function
        self.assertEqual(self.prune('chrOther', segments, links, {'chrUn_a', 'chrUn_b'})[1][:2], (segs, kept_links))

    def test_named_chromosome_never_reads_the_fasta(self):
        def bin_ref_contigs():
            raise AssertionError('the bin reference fasta is only for the --otherContig lump')
        self.assertEqual(self.prune('chr4', bin_ref_contigs=bin_ref_contigs)[1][0], self.KEPT_SEGMENTS)
        self.assertEqual(prune_own_contigs('chr4', ['chr4', 'chr10'], bin_ref_contigs), {'chr4'})

    def test_no_own_contig_is_an_error(self):
        with self.assertRaisesRegex(RuntimeError, 'nothing to prune it back to'):
            self.prune('chrOther', bin_ref_contigs={'chrUn_x'})
        with self.assertRaises(RuntimeError):
            self.prune('chr5')

    def test_nothing_to_prune(self):
        # only chr4's own material: nothing goes, and the job leaves the graphs as they were
        segments = [s for s in self.SEGMENTS if s[0] in ('s1', 's2', 's3', 's4', 's9')]
        links = [l for l in self.LINKS if l[0] in ('s1', 's2', 's9') and l[1] in ('s2', 's3', 's9')]
        scan = scan_rgfa(self.write_gfa(segments, links))
        prune = prune_other_contigs_plan(scan, {'chr4'})
        self.assertFalse(prune.dropped())
        self.assertEqual(prune_other_contigs_report(prune.summary)[-1], 'kept 5 of 5 segments (0 bp dropped) and 5 of 5 links')

    def test_not_an_rgfa(self):
        path = os.path.join(self.work_dir, 'plain.gfa')
        with open(path, 'w') as gfa_file:
            gfa_file.write('S\ts1\tACGT\nS\ts2\tACGT\nL\ts1\t+\ts2\t+\t0M\n')
        with self.assertRaisesRegex(RuntimeError, 'Segment s1 .* has no SR:i: tag: it is not an rGFA'):
            scan_rgfa(path)
        two_segments = rgfa_segment('s1', 'ACGT', 'id=CHM13|chr4', 0) + rgfa_segment('s2', 'ACGT', 'id=CHM13|chr4', 0)
        # a link without its rank could be a sample link onto another chromosome, which the prune
        # would then keep, so it is refused as a segment without one is
        with open(path, 'w') as gfa_file:
            gfa_file.write(two_segments + 'L\ts1\t+\ts2\t+\t0M\n')
        with self.assertRaisesRegex(RuntimeError, 'Link s1 \\+ s2 \\+ .* has no SR:i: tag: it is not an rGFA'):
            scan_rgfa(path)
        # a path can name a segment the prune drops, and the prune does not know how to cut one
        for path_line in ('W\tHG1\t1\tctg1\t0\t8\t>s1>s2\n', 'P\tHG1#1#ctg1\ts1+,s2+\t*\n'):
            with open(path, 'w') as gfa_file:
                gfa_file.write(two_segments + rgfa_link('s1', 's2', 0) + path_line)
            with self.assertRaisesRegex(RuntimeError, 'has {} lines'.format(path_line[0])):
                scan_rgfa(path)
        with open(path, 'w') as gfa_file:
            gfa_file.write(rgfa_segment('s1', 'ACGT', 'id=CHM13|chr4', 0) + rgfa_link('s1', 's2', 0))
        with self.assertRaisesRegex(RuntimeError, 'have no S line'):
            scan_rgfa(path)

    def test_config(self):
        root = ET.parse(os.path.join(cactusRootPath(), 'cactus_progressive_config.xml')).getroot()
        node = root.find('graphmap_split')
        # on in the default config, and in one without the attribute
        self.assertEqual(node.get('pruneOtherContigs'), '1')
        self.assertTrue(prune_other_contigs_enabled(root))
        del node.attrib['pruneOtherContigs']
        self.assertTrue(prune_other_contigs_enabled(root))
        for value in ('0', 'false', 'False'):
            node.set('pruneOtherContigs', value)
            self.assertFalse(prune_other_contigs_enabled(root))
        node.set('pruneOtherContigs', 'true')
        self.assertTrue(prune_other_contigs_enabled(root))
        # it applies to every chromosome or none: a list of them is refused rather than read as on
        for value in ('chr4,chr21', ''):
            node.set('pruneOtherContigs', value)
            with self.assertRaisesRegex(RuntimeError, 'must be 0 or 1'):
                prune_other_contigs_enabled(root)

    def test_memory(self):
        # 6x the longest reference contig, with a floor: CHM13 (chr1, 248 Mbp) is under the floor, a
        # genome with a 2 Gbp chromosome is not
        self.assertEqual(prune_other_contigs_memory(None), PRUNE_MIN_MEMORY)
        self.assertEqual(prune_other_contigs_memory(248387328), PRUNE_MIN_MEMORY)
        self.assertEqual(prune_other_contigs_memory(2 * 10**9), 12 * 10**9)


if __name__ == '__main__':
    unittest.main()
