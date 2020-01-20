#!/usr/bin/env python3

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
"""
"""

import unittest
import os
import sys

from sonLib.bioio import TestStatus
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger
from sonLib.bioio import system

from cactus.progressive.multiCactusTree import MultiCactusTree
from sonLib.nxnewick import NXNewick

class TestCase(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)
        self.mcTree1 = None
        self.mcTree1a = None
        self.mcTree2 = None
        self.__generateTrees()

    @TestStatus.shortLength
    def testSanity(self):
        parser = NXNewick()
        mcTree1 = MultiCactusTree(parser.parseString(self.tree1, addImpliedRoots = False))
        tree1String = NXNewick().writeString(mcTree1)
        self.assertEqual(tree1String, self.tree1)
        mcTree2 = MultiCactusTree(parser.parseString(self.tree2, addImpliedRoots = False))
        tree2String = NXNewick().writeString(mcTree2)
        self.assertEqual(tree2String, self.tree2)

    @TestStatus.shortLength
    def testSubtrees(self):
        roots1 = ["Anc0", "Anc1", "Anc2", "Anc3", "Anc4", "Anc5", "Anc6", "Anc7"]
        roots2 = ["Anc0", "Anc1", "Anc2", "Anc3", "Anc4", "Anc5"]

        subTree1_a3 = '(Anc7:0.025291,BABOON:0.044568)Anc3;'
        subTree2_a5 = '(monkey:100.8593,cat:47.14069)Anc5;'

        trueRoots = [roots1, roots2]
        trueSubtrees = [subTree1_a3, subTree2_a5]
        trees = [self.mcTree1, self.mcTree2]
        ancs = ["Anc3", "Anc5"]

        for tree, trueRoot, anc, trueSubtree in zip(trees, trueRoots, ancs, trueSubtrees):
            roots = tree.getSubtreeRootNames()
            self.assertEqual(sorted(roots), sorted(trueRoot))
            subtree = tree.extractSubTree(anc)
            subtree = NXNewick().writeString(subtree)
            self.assertEqual(subtree, trueSubtree)

    @TestStatus.shortLength
    def testAddSelf(self):
        trueSelf = '((((((((HUMAN:0.006969)HUMAN_self:0.006969,(CHIMP:0.009727)CHIMP_self:0.009727)Anc7:0.025291)Anc7_self:0.025291,(BABOON:0.044568)BABOON_self:0.044568)Anc3:0.11)Anc3_self:0.11,(((MOUSE:0.072818)MOUSE_self:0.072818,(RAT:0.081244)RAT_self:0.081244)Anc4:0.260342)Anc4_self:0.260342)Anc1:0.02326)Anc1_self:0.02326,(((((DOG:0.07)DOG_self:0.07,(CAT:0.07)CAT_self:0.07)Anc5:0.087381)Anc5_self:0.087381,(((PIG:0.06)PIG_self:0.06,(COW:0.06)COW_self:0.06)Anc6:0.104728)Anc6_self:0.104728)Anc2:0.04)Anc2_self:0.04)Anc0;'
        tree = MultiCactusTree(self.mcTree1)
        tree.nameUnlabeledInternalNodes()
        tree.computeSubtreeRoots()
        tree.addSelfEdges()
        treeString = NXNewick().writeString(tree)
        self.assertEqual(treeString, trueSelf)

    @TestStatus.shortLength
    def testAddOutgroup(self):
        trueOg = '((((HUMAN:0.006969,CHIMP:0.009727)Anc7:0.025291,BABOON:0.044568)Anc3:0.11,(MOUSE:0.072818,RAT:0.081244)Anc4:0.260342)Anc1:0.02326,((DOG:0.07,CAT:0.07)Anc5:0.087381,(PIG:0.06,COW:0.06)Anc6:0.104728)Anc2:0.04,outgroup:1.7)Anc0;'
        tree = MultiCactusTree(self.mcTree1)
        tree.nameUnlabeledInternalNodes()
        tree.computeSubtreeRoots()
        tree.addOutgroup("outgroup", 1.7)
        treeString = NXNewick().writeString(tree)
        self.assertEqual(treeString, trueOg)

        trueLeafOg = "(A:1.1,outgroup:1.1);"
        leafTreeString = "A;"
        parser = NXNewick()
        leafTree = MultiCactusTree(parser.parseString(leafTreeString, addImpliedRoots = False))
        leafTree.nameUnlabeledInternalNodes()
        leafTree.computeSubtreeRoots()
        leafTree.addOutgroup("outgroup", 2.2)
        leafTreeOutString = NXNewick().writeString(leafTree)
        self.assertEqual(leafTreeOutString, trueLeafOg)

    @TestStatus.shortLength
    def testExtractSpanningTree(self):
        """Tests whether extracting a binary spanning tree works correctly."""
        prevNewick1 = NXNewick().writeString(self.mcTree1)
        # Check a dead-simple spanning tree with 3 closely related leaves.
        spanHCB = self.mcTree1.extractSpanningTree(["HUMAN", "CHIMP", "BABOON"])
        # Check that the existing tree hasn't been modified (OK, a bit
        # silly, but just in case).
        self.assertEqual(NXNewick().writeString(self.mcTree1), prevNewick1)
        # Check the actual spanning tree.
        self.assertEqual(NXNewick().writeString(spanHCB), "((HUMAN:0.006969,CHIMP:0.009727)Anc7:0.025291,BABOON:0.044568)Anc3;")

        # Now test a more complicated tree, where we should remove as
        # many of the ancestors as possible (they will add extra
        # losses for no reason!).
        spanHCC = self.mcTree1.extractSpanningTree(["HUMAN", "CHIMP", "CAT"])
        self.assertEqual(NXNewick().writeString(self.mcTree1), prevNewick1)
        self.assertEqual(NXNewick().writeString(spanHCC), "((HUMAN:0.006969,CHIMP:0.009727)Anc7:0.158551,CAT:0.197381)Anc0;")

    @TestStatus.shortLength
    def testGetChildren(self):
        self.assertEqual(self.mcTree1.getChildNames('Anc6'), ['PIG', 'COW'])

    def __generateTrees(self):
        self.tree1 = '((((HUMAN:0.006969,CHIMP:0.009727):0.025291,BABOON:0.044568):0.11,(MOUSE:0.072818,RAT:0.081244):0.260342):0.02326,((DOG:0.07,CAT:0.07):0.087381,(PIG:0.06,COW:0.06):0.104728):0.04);'
        self.tree2 = '((raccoon:19.19959,bear:6.80041):0.846,((sea_lion:11.997,seal:12.003):7.52973,((monkey:100.8593,cat:47.14069):20.59201,weasel:18.87953):2.0946):3.87382,dog:25.46154);'
        parser = NXNewick()
        self.mcTree1 = MultiCactusTree(parser.parseString(self.tree1, addImpliedRoots = False))
        self.mcTree2 = MultiCactusTree(parser.parseString(self.tree2, addImpliedRoots = False))
        self.mcTree1.nameUnlabeledInternalNodes()
        self.mcTree2.nameUnlabeledInternalNodes()
        self.mcTree1.computeSubtreeRoots()
        self.mcTree2.computeSubtreeRoots()

def main():
    unittest.main()

if __name__ == '__main__':
    main()
