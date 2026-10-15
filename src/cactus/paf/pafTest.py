#!/usr/bin/env python3

"""
Copyright (C) 2009-2021 by Benedict Paten (benedictpaten@gmail.com)

Released under the MIT license, see LICENSE.txt
"""

import unittest
from sonLib.bioio import newickTreeParser
from cactus.paf.paf import *
from cactus.paf.last_scoring import apply_long_gap, apply_scores_to_config, abpoa_exact_depth, parse_train_file
from cactus.paf.last_scoring import write_pair_hmm_train_file, pick_train_partners
from cactus.paf.pair_hmm_training import maximise, combined_gap, LONG_GAP_OPEN, LONG_GAP_EXTEND
from toil.realtimeLogger import RealtimeLogger
import logging
import os
import tempfile
import xml.etree.ElementTree as ET

class TestCase(unittest.TestCase):
    def setUp(self):
        self.tree_string = "((simCow:0.1,simDog:0.2)anc1:0.3,(simChimp:0.4,simHuman:0.5)anc2:0.2)anc0:0.1;"
        self.tree = newickTreeParser(self.tree_string)

    def test_get_event_pairs(self):
        leaves = get_leaves(self.tree)
        self.assertEqual(6, len(list(get_event_pairs(self.tree, leaves))))
        distances = get_distances(self.tree)
        for species_a, species_b, distance in get_event_pairs(self.tree, leaves):
            self.assertTrue(species_a in leaves)
            self.assertTrue(species_b in leaves)
            self.assertEqual(distances[(species_a, species_b)], distance)

    def test_get_distances(self):
        distances = get_distances(self.tree)
        cow = get_node(self.tree, "simCow")
        human = get_node(self.tree, "simHuman")
        dog = get_node(self.tree, "simDog")
        anc2 = get_node(self.tree, "anc2")
        self.assertAlmostEqual(distances[(cow, human)], 0.1+0.3+0.5+0.2)
        self.assertAlmostEqual(distances[(cow, dog)], 0.1+0.2)
        self.assertAlmostEqual(distances[(dog, cow)], 0.1+0.2)
        self.assertAlmostEqual(distances[(anc2, human)], 0.5)
        self.assertAlmostEqual(distances[(anc2, cow)], 0.1+0.3+0.2)

    def test_get_node(self):
        for i in "simCow", "simDog", "simChimp", "simHuman", "anc1", "anc2", "anc0":
            self.assertEqual(i, get_node(self.tree, i).iD)

    def test_get_subtree_nodes(self):
        self.assertEqual({ "simCow", "simDog", "simChimp", "simHuman", "anc1", "anc2", "anc0" }, { i.iD for i in get_subtree_nodes(self.tree) })

    def test_get_leaves(self):
        self.assertEqual({ "simCow", "simDog", "simChimp", "simHuman" }, { i.iD for i in get_leaves(self.tree) })

class LastScoringTest(unittest.TestCase):
    """ how a last-train model is scaled for abPOA (cactus.paf.last_scoring) """
    ROWS = {'A': [5, -32, -26, -34], 'C': [-32, 7, -30, -26], 'G': [-26, -30, 7, -32], 'T': [-34, -26, -32, 5]}

    def model(self, gap_open, gap_extend):
        # the HPRC oct9 CHM13 chr2 model's matrix (that chromosome learned an extension of 2, the rest 1)
        d = {i: dict(zip('ACGT', self.ROWS[i])) for i in 'ACGT'}
        d['GAP-OPEN'], d['GAP-EXTEND'] = gap_open, gap_extend
        return d

    @staticmethod
    def gaps(d):
        return d['GAP-OPEN'], d['GAP-EXTEND'], d['GAP-OPEN-2'], d['GAP-EXTEND-2']

    @staticmethod
    def config():
        return ET.fromstring('<cactus><bar><poa partialOrderAlignmentWindow="10000" partialOrderAlignmentTrainedGapOpen2Factor="3" '
                             'partialOrderAlignmentTrainedGapExtension2Factor="30"/><minipoa/></bar></cactus>')

    @staticmethod
    def write_train(path, gap_open, gap_extend, rows):
        with open(path, 'w') as f:
            f.write('#last -a {0}\n#last -A {0}\n#last -b {1}\n#last -B {1}\n'.format(gap_open, gap_extend))
            f.write('       A      C      G      T\n')
            for i in 'ACGT':
                f.write('{} {}\n'.format(i, ' '.join(str(v) for v in rows[i])))

    def test_extension_one_unchanged(self):
        # the usual model: scaled 30x and nothing else, as before
        d = self.model(43, 1)
        self.assertEqual(apply_long_gap(d, 3, 30), 1)
        self.assertEqual(self.gaps(d), (1290, 30, 3870, 1))
        self.assertEqual(d['A']['C'], -960)

    def test_extension_two_reduced(self):
        # scaled 30x, every score is even: halved, the same model with abPOA's usual extension of 30
        d = self.model(43, 2)
        self.assertEqual(apply_long_gap(d, 3, 30), 2)
        self.assertEqual(self.gaps(d), (645, 30, 1935, 1))
        self.assertEqual([d['A'][j] for j in 'ACGT'], [75, -480, -390, -510])

    def test_reduction_is_exact(self):
        # whatever the extension, the result is the 30x model divided by a common factor
        for ext in range(1, 9):
            d, ref = self.model(43, ext), self.model(43, ext)
            div = apply_long_gap(d, 3, 30)
            self.assertEqual(div, {1: 1, 2: 2, 3: 3, 4: 2, 5: 5, 6: 6, 7: 1, 8: 2}[ext])
            self.assertEqual(self.gaps(d), (43 * 30 // div, ext * 30 // div, 43 * 90 // div, ext // div))
            for i in 'ACGT':
                for j in 'ACGT':
                    self.assertEqual(d[i][j] * div, ref[i][j] * 30)

    def test_config_gets_reduced_scores(self):
        root = self.config()
        apply_scores_to_config(self.model(43, 2), root, label='chr2')
        poa, minipoa = root.find('bar/poa'), root.find('bar/minipoa')
        self.assertEqual([poa.attrib['partialOrderAlignmentGapOpenPenalty1'], poa.attrib['partialOrderAlignmentGapExtensionPenalty1'],
                          poa.attrib['partialOrderAlignmentGapOpenPenalty2'], poa.attrib['partialOrderAlignmentGapExtensionPenalty2']],
                         ['645', '30', '1935', '1'])
        self.assertEqual(poa.attrib['partialOrderAlignmentSubMatrix'].split()[:5], ['75', '-480', '-390', '-510', '-390'])
        self.assertEqual((minipoa.attrib['minipoaGapOpenPenalty'], minipoa.attrib['minipoaGapExtensionPenalty']), ('645', '30'))

    def test_config_log(self):
        # dividing a model is exact, so nothing is warned about, and the override line keeps the prefix logs are grepped for
        RealtimeLogger.getLogger()  # outside a workflow its first use sets it to CRITICAL, which would undo assertLogs
        with self.assertLogs('toil-rt', level='INFO') as logs:
            apply_scores_to_config(self.model(43, 2), self.config(), label='chr2')
        self.assertEqual([r.getMessage() for r in logs.records if r.levelno >= logging.WARNING], [])
        override = [r.getMessage() for r in logs.records if r.getMessage().startswith(
            'Overriding abPOA scores with trained values: GapOpen 645; GapExtend 30; GapOpen2 1935; GapExtend2 1; SubMatrix 75 ')]
        self.assertEqual(len(override), 1)
        self.assertIn('chr2: learned gap extension 2; all scores multiplied by 30, then divided by their common factor 2', override[0])

    def test_exact_depth(self):
        # abPOA's own numbers for a trained model at the 10 kb window: far beyond any window's graph
        d = self.model(47, 1)
        apply_long_gap(d, 3, 30)
        mat = [d[i][j] for i in 'ACGT' for j in 'ACGT'] + [-990, 210, -990, 210, 150]
        self.assertGreater(abpoa_exact_depth(mat, *self.gaps(d), 10000), 100 * 10000)
        # four times the scores leaves less than 50 windows of depth, which is warned about
        self.assertLess(abpoa_exact_depth([4 * v for v in mat], *[4 * g for g in self.gaps(d)], 10000), 50 * 10000)
        # convex gaps whose two pieces have the same open + extension (110): the larger extension (60) still counts,
        # as in abPOA's simd_abpoa_score_range(), which gives this depth
        self.assertEqual(abpoa_exact_depth([10 if i % 6 == 0 else -20 for i in range(25)], 100, 10, 50, 60, 10000), 19509843)

    def test_unusable_model_rejected(self):
        # an extension below 1, a negative open, or any match score of 0 or less (here A and T, the others positive)
        zero_at = {i: [0 if i == j and i in 'AT' else v for j, v in zip('ACGT', self.ROWS[i])] for i in 'ACGT'}
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, 'model.train')
            self.write_train(path, 43, 1, self.ROWS)
            self.assertEqual(parse_train_file(path)['GAP-EXTEND'], 1)
            for gap_open, gap_extend, rows in [(43, 0, self.ROWS), (-1, 1, self.ROWS), (43, 1, zero_at)]:
                self.write_train(path, gap_open, gap_extend, rows)
                with self.assertRaises(RuntimeError):
                    parse_train_file(path)

class PairHmmTrainFileTest(unittest.TestCase):
    """ how a cPecan pair HMM becomes a .train file and abPOA's scores (cactus.paf.last_scoring) """

    @staticmethod
    def hmm(short_open=1e-4, short_extend=0.5, long_open=3e-5, long_extend=0.9565, identity=0.998):
        # five states, match and a short and a long gap in each sequence; transitions[from * 5 + to]
        T = [[0.0] * 5 for _ in range(5)]
        T[0] = [1 - 2 * short_open - 2 * long_open, short_open, short_open, long_open, long_open]
        for g in (1, 2):
            T[g][g], T[g][0] = short_extend, 1 - short_extend
        for g in (3, 4):
            T[g][g], T[g][0] = long_extend, 1 - long_extend
        # transitions twice as likely as each transversion, A/T a little commoner than C/G
        freq = {0: 0.27, 1: 0.23, 2: 0.23, 3: 0.27}
        match = [[freq[x] * identity if x == y else (1 - identity) * (2 if (x + y) % 2 == 0 else 1) / 16
                  for y in range(4)] for x in range(4)]
        total = sum(map(sum, match))
        emissions = [match[x][y] / total for x in range(4) for y in range(4)] + [1 / 16] * 64
        return [v for row in T for v in row], emissions

    def write_and_parse(self, tmp, transitions, emissions, **kwargs):
        path = os.path.join(tmp, 'model.train')
        write_pair_hmm_train_file(path, transitions, emissions, 'test', 3, 30, **kwargs)
        return parse_train_file(path)

    def test_trained_long_gap_is_second_piece(self):
        # a long state of 23bp on average is a real long gap: dearer to open, cheaper to extend, and kept as it is
        with tempfile.TemporaryDirectory() as tmp:
            d = self.write_and_parse(tmp, *self.hmm())
        self.assertGreater(d['GAP-OPEN-2'], d['GAP-OPEN'])
        self.assertLess(d['GAP-EXTEND-2'], d['GAP-EXTEND'])
        self.assertNotEqual(d['GAP-OPEN-2'], 3 * d['GAP-OPEN'])
        for i in 'ACGT':
            self.assertGreater(d[i][i], 0)
            for j in 'ACGT':
                self.assertEqual(d[i][j], d[j][i])
        # apply_scores_to_config leaves the two pieces as they are (but for a common factor), where a one-piece model's
        # second piece is made from its first
        root = LastScoringTest.config()
        pieces = LastScoringTest.gaps(dict(d))
        apply_scores_to_config(d, root)
        poa = root.find('bar/poa')
        got = [int(poa.attrib[k]) for k in ('partialOrderAlignmentGapOpenPenalty1', 'partialOrderAlignmentGapExtensionPenalty1',
                                            'partialOrderAlignmentGapOpenPenalty2', 'partialOrderAlignmentGapExtensionPenalty2')]
        factor = pieces[0] // got[0]
        self.assertEqual([g * factor for g in got], list(pieces))
        self.assertEqual(root.find('bar/minipoa').attrib['minipoaGapOpenPenalty2'], '0')

    def test_short_long_state_falls_back(self):
        # a "long" state of under 10bp on average is a second short gap: the first piece is both states as one, and
        # the second is made from it with the factors, as apply_long_gap makes one
        with tempfile.TemporaryDirectory() as tmp:
            d = self.write_and_parse(tmp, *self.hmm(long_extend=0.79))
            held = self.write_and_parse(tmp, *self.hmm(), long_gap_trained=False)
        self.assertEqual((d['GAP-OPEN-2'], d['GAP-EXTEND-2']), (3 * d['GAP-OPEN'], max(1, d['GAP-EXTEND'] // 30)))
        # held at cPecan's default, the long state is not trained, so the second piece is made the same way, from the
        # short state alone
        self.assertEqual((held['GAP-OPEN-2'], held['GAP-EXTEND-2']), (3 * held['GAP-OPEN'], max(1, held['GAP-EXTEND'] // 30)))
        self.assertGreater(held['GAP-EXTEND'], d['GAP-EXTEND'])

    def test_combined_gap_keeps_mean_length(self):
        transitions, _ = self.hmm(short_open=1e-4, short_extend=0.5, long_open=1e-4, long_extend=0.75)
        open_prob, extend_prob = combined_gap(transitions)
        self.assertAlmostEqual(open_prob, 2e-4)
        # gaps of 2 and 4bp on average, opened as often: 3bp
        self.assertAlmostEqual(1 / (1 - extend_prob), 3.0)

    def test_maximise_holds_long_gap(self):
        # expected counts with a long state that looks like a short one: held at cPecan's unless train_long_gap
        transitions, emissions = self.hmm(long_extend=0.5)
        held, _ = maximise([v * 1000 for v in transitions], [v * 1000 for v in emissions])
        self.assertAlmostEqual(held[0 * 5 + 3], LONG_GAP_OPEN)
        self.assertAlmostEqual(held[3 * 5 + 3], LONG_GAP_EXTEND)
        trained, _ = maximise([v * 1000 for v in transitions], [v * 1000 for v in emissions], train_long_gap=True)
        self.assertAlmostEqual(trained[3 * 5 + 3], 0.5, places=3)
        for t in (held, trained):
            for i in range(5):
                self.assertAlmostEqual(sum(t[i * 5:i * 5 + 3]) + sum(t[i * 5 + 3:i * 5 + 5]), 1.0 if i < 3 else
                                       1.0 + (LONG_GAP_OPEN if t is held else 0.0), places=6)

    def test_pick_train_partners(self):
        # the candidates are the genomes of typical size, in seq_order (by mash distance), and the partners the ones
        # nearest the middle of them; a fragment and an oversized assembly are never picked
        class Fasta:
            def __init__(self, size):
                self.size = size
        sizes = {'ref': 10**8, 'a': 10**8, 'b': 10**8, 'tiny': 10**5, 'c': 10**8, 'd': 10**8, 'huge': 5 * 10**8,
                 'e': 10**8, 'f': 10**8}
        seq_id_map = {name: Fasta(size) for name, size in sizes.items()}
        order = ['ref', 'a', 'b', 'tiny', 'c', 'd', 'huge', 'e', 'f']
        self.assertEqual(pick_train_partners('ref', order, seq_id_map, 4), ['b', 'c', 'd', 'e'])
        self.assertEqual(pick_train_partners('ref', order, seq_id_map, 1), ['d'])
        self.assertEqual(pick_train_partners('ref', order, seq_id_map, 10), ['a', 'b', 'c', 'd', 'e', 'f'])
        self.assertEqual(pick_train_partners('ref', ['ref', 'tiny'], seq_id_map, 4), [])

if __name__ == '__main__':
    unittest.main()