#!/usr/bin/env python3

"""
Unit tests for cactus-phast's phyloFit underflow screen (fit_neutral_model): the SS-file helpers,
and the pilot / screen / check / refit logic with phyloFit replaced by a stand-in that says which
columns underflow.  Fast and offline -- no phast binaries.

The end-to-end workflow is smoke-tested in test/evolverTest.py.
"""

import argparse
import os
import shutil
import tempfile
import unittest
from collections import Counter
from unittest import mock

import cactus.maf.cactus_phast as cactus_phast
from cactus.maf.cactus_phast import (
    read_ss_header, read_ss_counts, iter_ss_tuples, write_ss_subset, sample_ss_counts,
    fit_neutral_model, column_log2_probs, read_colprobs, below)


def write_ss(path, names, columns):
    """ Write an unordered, single-column SS file (as msa_view --aggregate --unordered-ss does)
    from a list of column strings. """
    tuples = sorted(Counter(columns).items())
    with open(path, 'w') as f:
        f.write('NSEQS = {}\nLENGTH = {}\nTUPLE_SIZE = 1\nNTUPLES = {}\nNAMES = {}\n'
                'ALPHABET = ACGT\nNCATS = -1\n\n'.format(len(names), len(columns), len(tuples),
                                                       ','.join(names)))
        for i, (col, n) in enumerate(tuples):
            f.write('{}\t{}\t{}\n'.format(i, col, n))


def read_columns(path):
    """ The columns of an SS file with their counts. """
    return Counter({col: n for _, col, n in iter_ss_tuples(path)})


class TestSSFiles(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.ss = os.path.join(self.dir, 'in.ss')
        write_ss(self.ss, ['a', 'b', 'c'], ['AAA'] * 5 + ['ACG'] * 3 + ['TTT'] * 2 + ['GCA'])

    def tearDown(self):
        shutil.rmtree(self.dir)

    def test_counts(self):
        # tuples are sorted by column: AAA, ACG, GCA, TTT
        self.assertEqual(read_ss_counts(self.ss), [5, 3, 1, 2])
        lines, fields = read_ss_header(self.ss)
        self.assertEqual(fields['LENGTH'], '11')
        self.assertEqual(lines[-1], '')

    def test_subset_drops_and_renumbers(self):
        out = os.path.join(self.dir, 'out.ss')
        self.assertEqual(write_ss_subset(self.ss, out, [5, 0, 1, 2]), (8, 3))
        self.assertEqual(read_ss_counts(out), [5, 1, 2])
        self.assertEqual(read_columns(out), Counter({'AAA': 5, 'GCA': 1, 'TTT': 2}))
        _, fields = read_ss_header(out)
        self.assertEqual((fields['LENGTH'], fields['NTUPLES'], fields['NAMES']), ('8', '3', 'a,b,c'))

    def test_sample(self):
        counts = read_ss_counts(self.ss)
        s = sample_ss_counts(counts, 6, 1)
        self.assertEqual(sum(s), 6)
        self.assertTrue(all(0 <= x <= c for x, c in zip(s, counts)))
        self.assertEqual(s, sample_ss_counts(counts, 6, 1))
        self.assertEqual(sample_ss_counts(counts, 11, 1), counts)

    def test_bad_counts(self):
        bad = os.path.join(self.dir, 'bad.ss')
        with open(self.ss) as f, open(bad, 'w') as out:
            out.write(f.read().replace('LENGTH = 11', 'LENGTH = 12'))
        with self.assertRaises(RuntimeError):
            read_ss_counts(bad)


class FakePhyloFit:
    """ Stands in for run_phyloFit and column_log2_probs.  A "model" is a file naming the SS it
    was fit to; under any model, the columns in `junk` have a log2 likelihood of -5000 (as an
    underflowed one would) and the rest -50.  With always_bad, tuple 0 of whatever is scored is
    also -5000.  Records the SS of every fit. """

    def __init__(self, junk, always_bad=False):
        self.junk = set(junk)
        self.always_bad = always_bad
        self.fits = []

    def run_phyloFit(self, options, ss_path, tree_path, out_root, label='phyloFit'):
        self.fits.append((label, ss_path))
        with open(out_root + '.mod', 'w') as f:
            f.write(ss_path + '\n')
        return out_root + '.mod'

    def column_log2_probs(self, options, ss_path, mod_path, out_root):
        probs = [-5000.0 if col in self.junk else -50.0 for _, col, _ in iter_ss_tuples(ss_path)]
        if self.always_bad:
            probs[0] = -5000.0
        return probs


def all_columns(n_species):
    """ Every column of n_species bases over ACGT, once each. """
    cols = ['']
    for _ in range(n_species):
        cols = [c + b for c in cols for b in 'ACGT']
    return cols


class TestUnderflowScreen(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.work = os.path.join(self.dir, 'work')
        os.mkdir(self.work)
        names = ['a', 'b', 'c', 'd', 'e']
        # 1024 distinct columns, 2048 in all, plus 2 junk ones (0.1%)
        self.clean_cols = all_columns(5) * 2
        self.junk = ['ACGTN', 'TGCAN']
        self.ss = os.path.join(self.dir, 'all.4d.ss')
        write_ss(self.ss, names, self.clean_cols + self.junk)
        self.clean_ss = os.path.join(self.dir, 'clean.4d.ss')
        write_ss(self.clean_ss, names, self.clean_cols)

    def tearDown(self):
        shutil.rmtree(self.dir)

    def fit(self, ss, fake, keep=False, pilot_columns=1000):
        opts = argparse.Namespace(keepUnderflowColumns=keep, phyloFitCores=1)
        with mock.patch.object(cactus_phast, 'run_phyloFit', fake.run_phyloFit), \
             mock.patch.object(cactus_phast, 'column_log2_probs', fake.column_log2_probs), \
             mock.patch.object(cactus_phast, 'PILOT_COLUMNS', pilot_columns):
            mod = fit_neutral_model(self.work, opts, ss, 'tree.nh')
        with open(mod) as f:
            return f.read().strip()  # the SS the returned model was fit to

    def assertClean(self, ss_path):
        self.assertEqual(read_columns(ss_path), Counter(self.clean_cols))

    def test_nothing_to_screen(self):
        fake = FakePhyloFit(self.junk)
        self.assertEqual(self.fit(self.clean_ss, fake), self.clean_ss)
        self.assertEqual(fake.fits, [('phyloFit', self.clean_ss)])

    def test_keep_underflow_columns(self):
        fake = FakePhyloFit(self.junk)
        self.assertEqual(self.fit(self.ss, fake, keep=True), self.ss)
        self.assertEqual(len(fake.fits), 1)

    def test_check_refits_without_underflowing_columns(self):
        # 2050 columns is under PILOT_MIN_FACTOR pilots: fit everything, the check catches the
        # junk, refit
        fake = FakePhyloFit(self.junk)
        fit_ss = self.fit(self.ss, fake)
        self.assertEqual([label for label, _ in fake.fits], ['phyloFit', 'phyloFit'])
        self.assertEqual(fake.fits[0][1], self.ss)
        self.assertEqual(fit_ss, fake.fits[1][1])
        self.assertClean(fit_ss)

    def test_pilot_screens_underflowing_columns(self):
        fake = FakePhyloFit(self.junk)
        fit_ss = self.fit(self.ss, fake, pilot_columns=100)
        pilot_fits = [ss for label, ss in fake.fits if label == 'phyloFit-pilot']
        self.assertTrue(pilot_fits)
        self.assertEqual(read_ss_header(pilot_fits[0])[1]['LENGTH'], '100')
        # the screen left the junk out, so the real fit ran once, on the clean columns
        self.assertEqual([ss for label, ss in fake.fits if label == 'phyloFit'], [fit_ss])
        self.assertEqual(os.path.basename(fit_ss), 'train.ss')
        self.assertClean(fit_ss)

    def test_too_much_left_out(self):
        # junk that's 25% of the columns is not a few extreme columns
        fake = FakePhyloFit(all_columns(5)[:256])
        with self.assertRaises(RuntimeError):
            self.fit(self.clean_ss, fake)
        with self.assertRaises(RuntimeError):
            self.fit(self.clean_ss, fake, pilot_columns=100)

    def test_gives_up(self):
        # a different tuple underflows every time: refit MAX_UNDERFLOW_REFITS times, then fail
        fake = FakePhyloFit([], always_bad=True)
        with self.assertRaises(RuntimeError):
            self.fit(self.clean_ss, fake)
        self.assertEqual(len(fake.fits), cactus_phast.MAX_UNDERFLOW_REFITS + 1)
        # only the input and the latest refit's subset are left
        self.assertEqual(sorted(f for f in os.listdir(self.work) if f.endswith('.ss')),
                         ['phyloFit.{}.ss'.format(cactus_phast.MAX_UNDERFLOW_REFITS)])


class TestColumnLog2Probs(unittest.TestCase):
    """ column_log2_probs with phyloFit replaced: it must score a unit-count copy of the SS
    (LENGTH == NTUPLES; otherwise phyloFit reads past its buffer) and map line j to tuple j. """

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.ss = os.path.join(self.dir, 'in.ss')
        write_ss(self.ss, ['a', 'b', 'c'], ['AAA'] * 5000 + ['ACG'] * 3 + ['TTT'] * 2 + ['GCA'])

    def tearDown(self):
        shutil.rmtree(self.dir)

    def fake_phyloFit(self, parameters, **kwargs):
        scored = parameters[-1]
        out_root = parameters[parameters.index('--out-root') + 1]
        _, fields = read_ss_header(scored)
        self.scored_fields = fields
        self.scored_columns = read_columns(scored)
        with open(out_root + '.colprobs', 'w') as f:
            for j, (_, col, _) in enumerate(iter_ss_tuples(scored)):
                f.write('{}\t{:.6f}\n'.format(j, -float(len(set(col)))))

    def test_scores_unit_count_copy(self):
        with mock.patch.object(cactus_phast, 'cactus_call', self.fake_phyloFit):
            probs = column_log2_probs(argparse.Namespace(phyloFitCores=1), self.ss, 'm.mod',
                                      os.path.join(self.dir, 'out'))
        self.assertEqual(self.scored_fields['LENGTH'], self.scored_fields['NTUPLES'])
        self.assertEqual(self.scored_columns, Counter({'AAA': 1, 'ACG': 1, 'GCA': 1, 'TTT': 1}))
        # tuples are AAA, ACG, GCA, TTT, scored -(distinct bases)
        self.assertEqual(probs, [-1.0, -3.0, -3.0, -1.0])
        self.assertFalse(os.path.exists(os.path.join(self.dir, 'out.unit.ss')))

    def test_read_colprobs_checks_lines(self):
        p = os.path.join(self.dir, 'x.colprobs')
        with open(p, 'w') as f:
            f.write('0\t-1.5\n1\t-999999999.000000\n2\tnan\n')
        probs = read_colprobs(p, 3)
        self.assertEqual(probs[:2], [-1.5, -999999999.0])
        self.assertTrue(below(probs[2], -700))
        self.assertTrue(below(probs[1], -700))
        self.assertFalse(below(probs[0], -700))
        with self.assertRaises(RuntimeError):
            read_colprobs(p, 4)
        with open(p, 'w') as f:
            f.write('0\t-1.5\n2\t-1.5\n')
        with self.assertRaises(RuntimeError):
            read_colprobs(p, 2)


if __name__ == '__main__':
    unittest.main()
