#!/usr/bin/env python3

"""
Copyright (C) 2009-2021 by Benedict Paten (benedictpaten@gmail.com)

Released under the MIT license, see LICENSE.txt
"""

import unittest
from sonLib.bioio import newickTreeParser
from paf import *

class TestCase(unittest.TestCase):
    def setUp(self):
        self.tree_string = "((simCow:0.1,simDog:0.2)anc1:0.3,(simChimp:0.4,simHuman:0.5)anc2:0.2)anc0:0.1;"
        self.tree = newickTreeParser(self.tree_string)

    def test_get_leaf_event_pairs(self):
        self.assertEqual(6, len(list(get_leaf_event_pairs(self.tree))))
        distances = get_distances(self.tree)
        leaves = get_leaves(self.tree)
        for species_a, species_b, distance in get_leaf_event_pairs(self.tree):
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

if __name__ == '__main__':
    unittest.main()