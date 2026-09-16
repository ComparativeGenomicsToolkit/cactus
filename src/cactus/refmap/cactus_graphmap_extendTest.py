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
"""

import gzip
import os
import re
import tempfile
import unittest

from cactus.refmap.cactus_graphmap import (
    gaf_from_pansn, gaf_to_pansn, pansn_to_event_map, split_gaf_file_by_event, trim_unstable_gaf)


# a stable GAF line as minigraph writes it: query name, then a path of stable segments with
# orientation marks, then the tags.  only columns 1 and 6 carry sequence names
def gaf_line(query, path, tags='60\ttp:A:P\tcm:i:100'):
    return '{}\t1000\t0\t900\t+\t{}\t2000\t10\t910\t880\t900\t{}\n'.format(query, path, tags)


class TestPansnToEventMap(unittest.TestCase):

    def test_haploid_and_diploid(self):
        self.assertEqual(pansn_to_event_map({'S288C', 'HG002.1', 'HG002.2'}),
                         {'S288C#0': 'S288C', 'HG002#1': 'HG002.1', 'HG002#2': 'HG002.2'})

    def test_explicit_hap_zero_beats_bare_name(self):
        # event_to_pansn_prefix maps both S288C and S288C.0 onto S288C#0.  whichever wins, the map
        # has to be single valued or the round trip would be ambiguous
        prefix_map = pansn_to_event_map({'S288C', 'S288C.0'})
        self.assertEqual(len(prefix_map), 1)
        self.assertIn(prefix_map['S288C#0'], ('S288C', 'S288C.0'))

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
        names = {'HG002.1', 'HG002.2', 'GRCh38'}
        gaf = gaf_line('id=HG002.2|chr1', '<id=GRCh38|chr1>id=HG002.1|chr1<id=GRCh38|chr1')
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
        names = {'HG002.1', 'GRCh38'}
        gaf = gaf_line('id=HG002.1|chr1#0', '>id=GRCh38|chr1#0')
        pansn, back = self.round_trip(names, gaf)
        self.assertIn('HG002#1#chr1#0', pansn)
        self.assertEqual(back, gaf)

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
        names = {'HG002.1', 'HG002.2'}
        h1 = gaf_line('HG002#1#chr1', '>GRCh38#0#chr1')
        h2 = gaf_line('HG002#2#chr1', '>GRCh38#0#chr1')
        contents, _ = self.split(h1 + h2, names)
        self.assertEqual(contents, {'HG002.1': h1, 'HG002.2': h2})

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


if __name__ == '__main__':
    unittest.main()
