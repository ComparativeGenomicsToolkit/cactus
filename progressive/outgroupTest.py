#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
"""
"""

import unittest
import os
import sys
import copy
import xml.etree.ElementTree as ET
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger
from sonLib.bioio import system

from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.progressive.outgroup import GreedyOutgroup

from sonLib.nxnewick import NXNewick
from sonLib.nxtreeTest import randomTreeSet

class TestCase(unittest.TestCase):
    
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.trees = randomTreeSet()
        self.mcTrees = []
        for tree in self.trees:
            if tree.size() < 500:
                mcTree = MultiCactusTree(tree, tree.degree())
                for i in mcTree.breadthFirstTraversal():
                    mcTree.setName(i, "Node%s" % str(i))
                mcTree.computeSubtreeRoots()
                self.mcTrees.append(mcTree)
    
    def testJustLeaves(self):
        tree = '((((HUMAN:0.006969,CHIMP:0.009727)Anc7:0.025291,BABOON:0.044568)Anc6:0.11,(MOUSE:0.072818,RAT:0.081244)Anc5:0.260342)Anc4:0.023260,((DOG:0.07,CAT:0.07)Anc3:0.087381,(PIG:0.06,COW:0.06)Anc2:0.104728)Anc1:0.04)Anc0;'
        mcTree = MultiCactusTree(NXNewick().parseString(tree))
        mcTree.computeSubtreeRoots()        
        og = GreedyOutgroup()
        og.importTree(mcTree)
        candidates = set([mcTree.getName(x) for x in mcTree.getLeaves()])
        og.greedy(candidateSet=candidates, candidateChildFrac=2.)
        assert og.ogMap['Anc1'][0] == 'HUMAN'
        assert og.ogMap['Anc2'][0] in ['CAT', 'DOG']
        assert og.ogMap['Anc3'][0] in ['PIG', 'COW']
        assert og.ogMap['Anc4'][0] in ['CAT', 'DOG']
        assert og.ogMap['Anc5'][0] == 'HUMAN'
        assert og.ogMap['Anc6'][0] in ['CAT', 'DOG']
        assert og.ogMap['Anc7'][0] == 'BABOON'

    def testCandidates(self):
        tree = '((((HUMAN:0.006969,CHIMP:0.009727)Anc7:0.025291,BABOON:0.044568)Anc6:0.11,(MOUSE:0.072818,RAT:0.081244)Anc5:0.260342)Anc4:0.023260,((DOG:0.07,CAT:0.07)Anc3:0.087381,(PIG:0.06,COW:0.06)Anc2:0.104728)Anc1:0.04)Anc0;'
        mcTree = MultiCactusTree(NXNewick().parseString(tree))
        mcTree.computeSubtreeRoots()        
        og = GreedyOutgroup()
        og.importTree(mcTree)
        candidates = set(['HUMAN', 'CHIMP', 'RAT'])
        og.greedy(candidateSet=candidates, candidateChildFrac=0.5)
        assert og.ogMap['Anc1'][0] == 'Anc4'
        assert og.ogMap['Anc2'][0] == 'Anc4'
        assert og.ogMap['Anc3'][0] == 'Anc4'
        assert 'Anc4' not in og.ogMap
        assert og.ogMap['Anc5'][0] in ['HUMAN', 'CHIMP', 'Anc6', 'Anc7']
        assert og.ogMap['Anc6'][0] in ['Anc5', 'MOUSE', 'RAT']
        assert og.ogMap['Anc7'][0] in ['Anc5', 'MOUSE', 'RAT']

        og = GreedyOutgroup()
        og.importTree(mcTree)
        candidates = set(['HUMAN', 'CHIMP', 'RAT'])
        candidateFrac = 1        
        og.greedy(candidateSet=candidates, candidateChildFrac=1.0)
        assert og.ogMap['Anc1'][0] == 'Anc7'
        assert og.ogMap['Anc2'][0] == 'Anc7'
        assert og.ogMap['Anc3'][0] == 'Anc7'
        assert 'Anc4' not in og.ogMap
        assert og.ogMap['Anc5'][0] in ['HUMAN', 'CHIMP', 'Anc7']
        assert og.ogMap['Anc6'][0] == 'RAT'
        assert og.ogMap['Anc7'][0] == 'RAT'

    
            
    def testGeneralBetterThanLeaves(self):
        for tree in self.mcTrees:
            og1 = GreedyOutgroup()
            og1.importTree(tree)
            candidates = set([tree.getName(x) for x in tree.getLeaves()])
            og1.greedy(candidateSet=candidates, candidateChildFrac=2.)
            og2 = GreedyOutgroup()
            og2.importTree(tree)
            og2.greedy(candidateSet=None)
            
            for i in og1.ogMap:
                assert i in og2.ogMap
                dist1 = og1.ogMap[i][1]
                dist2 = og2.ogMap[i][1]
                assert dist2 <= dist1
                
    def testGeneralConstrainedBetterThanLeaves(self):
        for tree in self.mcTrees:
            og1 = GreedyOutgroup()
            og1.importTree(tree)
            candidates = set([tree.getName(x) for x in tree.getLeaves()])
            og1.greedy(candidateSet=candidates, candidateChildFrac=2.)
            og2 = GreedyOutgroup()
            og2.importTree(tree)
            og2.greedy(candidateSet=None, threshold=2)
            
            for i in og1.ogMap:
                assert i in og2.ogMap
                dist1 = og1.ogMap[i][1]
                dist2 = og2.ogMap[i][1]
                assert dist2 <= dist1
    
def main():
    unittest.main()
    
if __name__ == '__main__':
    main()
