#!/usr/bin/env python3

"""
Copyright (C) 2009-2021 by Benedict Paten (benedictpaten@gmail.com)

Released under the MIT license, see LICENSE.txt
"""

import unittest
from sonLib.bioio import newickTreeParser
from cactus.paf.paf import *
from cactus.paf.last_scoring import apply_long_gap, apply_scores_to_config, abpoa_exact_depth, parse_train_file
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

if __name__ == '__main__':
    unittest.main()