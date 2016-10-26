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
import random
from operator import itemgetter
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger
from sonLib.bioio import system

from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.progressive.outgroup import GreedyOutgroup, DynamicOutgroup

from sonLib.nxnewick import NXNewick
from sonLib.nxtreeTest import randomTreeSet

class TestCase(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)
        self.trees = randomTreeSet()
        self.mcTrees = []
        self.tempDir = getTempDirectory(os.getcwd())
        self.tempFa = os.path.join(self.tempDir, "seq.fa")
        with open(self.tempFa, "w") as f:
            f.write(">temp\nNNNNNNNCNNNNAAAAAAAAAAAAAAANNNNNNN\n")
        self.dummySeqMaps = []
        for tree in self.trees:
            if tree.size() < 500:
                mcTree = MultiCactusTree(tree, tree.degree())
                seqMap = dict()
                for i in mcTree.breadthFirstTraversal():
                    mcTree.setName(i, "Node%s" % str(i))
                    seqMap["Node%s" % str(i)] = self.tempFa
                mcTree.computeSubtreeRoots()
                self.mcTrees.append(mcTree)
                self.dummySeqMaps.append(seqMap)
                
        seqLens = dict()
        seqLens["HUMAN"] = 57553
        seqLens["CHIMP"] = 57344
        seqLens["BABOON"] = 58960
        seqLens["MOUSE"] = 32750
        seqLens["RAT"] = 38436
        seqLens["DOG"] = 54187
        seqLens["CAT"] = 50283
        seqLens["PIG"] = 54843
        seqLens["COW"] = 55508
        self.blanchetteSeqMap = dict()
        for event, seqLen in seqLens.items():
            p = os.path.join(self.tempDir, event +".fa")
            with open(p, "w") as f:
                f.write(">%s\n" % event)
                f.write(''.join(['A'] * seqLen))
                f.write('\n')
            self.blanchetteSeqMap[event] = p

    def tearDown(self):
        unittest.TestCase.tearDown(self)
        system("rm -rf %s" % self.tempDir)

    def testJustLeaves(self):
        tree = '((((HUMAN:0.006969,CHIMP:0.009727)Anc7:0.025291,BABOON:0.044568)Anc6:0.11,(MOUSE:0.072818,RAT:0.081244)Anc5:0.260342)Anc4:0.023260,((DOG:0.07,CAT:0.07)Anc3:0.087381,(PIG:0.06,COW:0.06)Anc2:0.104728)Anc1:0.04)Anc0;'
        mcTree = MultiCactusTree(NXNewick().parseString(tree, addImpliedRoots = False))
        mcTree.computeSubtreeRoots()
        og = GreedyOutgroup()
        og.importTree(mcTree)
        candidates = set([mcTree.getName(x) for x in mcTree.getLeaves()])
        og.greedy(candidateSet=candidates, candidateChildFrac=2.)
        assert og.ogMap['Anc1'][0][0] == 'HUMAN'
        assert og.ogMap['Anc2'][0][0] in ['CAT', 'DOG']
        assert og.ogMap['Anc3'][0][0] in ['PIG', 'COW']
        assert og.ogMap['Anc4'][0][0] in ['CAT', 'DOG']
        assert og.ogMap['Anc5'][0][0] == 'HUMAN'
        assert og.ogMap['Anc6'][0][0] in ['CAT', 'DOG']
        assert og.ogMap['Anc7'][0][0] == 'BABOON'

    def testCandidates(self):
        tree = '((((HUMAN:0.006969,CHIMP:0.009727)Anc7:0.025291,BABOON:0.044568)Anc6:0.11,(MOUSE:0.072818,RAT:0.081244)Anc5:0.260342)Anc4:0.023260,((DOG:0.07,CAT:0.07)Anc3:0.087381,(PIG:0.06,COW:0.06)Anc2:0.104728)Anc1:0.04)Anc0;'
        mcTree = MultiCactusTree(NXNewick().parseString(tree, addImpliedRoots = False))
        mcTree.computeSubtreeRoots()
        og = GreedyOutgroup()
        og.importTree(mcTree)
        candidates = set(['HUMAN', 'CHIMP', 'RAT'])
        og.greedy(candidateSet=candidates, candidateChildFrac=0.5)
        assert og.ogMap['Anc1'][0][0] == 'Anc4'
        assert og.ogMap['Anc2'][0][0] == 'Anc4'
        assert og.ogMap['Anc3'][0][0] == 'Anc4'
        assert 'Anc4' not in og.ogMap
        assert og.ogMap['Anc5'][0][0] in ['HUMAN', 'CHIMP', 'Anc6', 'Anc7']
        assert og.ogMap['Anc6'][0][0] in ['Anc5', 'MOUSE', 'RAT']
        assert og.ogMap['Anc7'][0][0] in ['Anc5', 'MOUSE', 'RAT']

        og = GreedyOutgroup()
        og.importTree(mcTree)
        candidates = set(['HUMAN', 'CHIMP', 'RAT'])
        candidateFrac = 1
        og.greedy(candidateSet=candidates, candidateChildFrac=1.0)
        assert og.ogMap['Anc1'][0][0] == 'Anc7'
        assert og.ogMap['Anc2'][0][0] == 'Anc7'
        assert og.ogMap['Anc3'][0][0] == 'Anc7'
        assert 'Anc4' not in og.ogMap
        assert og.ogMap['Anc5'][0][0] in ['HUMAN', 'CHIMP', 'Anc7']
        assert og.ogMap['Anc6'][0][0] == 'RAT'
        assert og.ogMap['Anc7'][0][0] == 'RAT'

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
                dist1 = og1.ogMap[i][0][1]
                dist2 = og2.ogMap[i][0][1]
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
                dist1 = og1.ogMap[i][0][1]
                dist2 = og2.ogMap[i][0][1]
                assert dist2 <= dist1

    def testMultipleOutgroups(self):
        tree = '((((HUMAN:0.006969,CHIMP:0.009727)Anc7:0.025291,BABOON:0.044568)Anc6:0.11,(MOUSE:0.072818,RAT:0.081244)Anc5:0.260342)Anc4:0.023260,((DOG:0.07,CAT:0.07)Anc3:0.087381,(PIG:0.06,COW:0.06)Anc2:0.104728)Anc1:0.04)Anc0;'
        mcTree = MultiCactusTree(NXNewick().parseString(tree, addImpliedRoots = False))
        mcTree.computeSubtreeRoots()
        og = GreedyOutgroup()
        og.importTree(mcTree)
        og.greedy(candidateChildFrac=0.5, maxNumOutgroups=3)
        # make sure all entries have <= 3 outgroups.
        assert all(map(lambda x: len(x) <= 3, og.ogMap.values()))
        # and for all entries, the closest must be first.
        assert all(map(lambda x: x == sorted(x, key=itemgetter(1)),
                       og.ogMap.values()))
        # ordering is important!
        assert map(itemgetter(0), og.ogMap['Anc4']) == ['Anc1']
        assert map(itemgetter(0), og.ogMap['Anc7']) == ['BABOON', 'Anc1',
                                                        'Anc5']
        # We avoid cycles, and choose post-order first, so this only
        # uses leaves.
        assert map(itemgetter(0), og.ogMap['Anc1']) == ['HUMAN', 'CHIMP',
                                                        'BABOON']

    def testMultipleOutgroupsJustLeaves(self):
        tree = '((((HUMAN:0.006969,CHIMP:0.009727)Anc7:0.025291,BABOON:0.044568)Anc6:0.11,(MOUSE:0.072818,RAT:0.081244)Anc5:0.260342)Anc4:0.023260,((DOG:0.07,CAT:0.07)Anc3:0.087381,(PIG:0.06,COW:0.06)Anc2:0.104728)Anc1:0.04)Anc0;'
        mcTree = MultiCactusTree(NXNewick().parseString(tree, addImpliedRoots = False))
        mcTree.computeSubtreeRoots()
        og = GreedyOutgroup()
        og.importTree(mcTree)
        candidates = set([mcTree.getName(x) for x in mcTree.getLeaves()])
        og.greedy(candidateSet=candidates, candidateChildFrac=2.,
                  maxNumOutgroups=3)
        # make sure all entries have <= 3 outgroups.
        assert all(map(lambda x: len(x) <= 3, og.ogMap.values()))
        # and for all entries, the closest must be first.
        assert all(map(lambda x: x == sorted(x, key=itemgetter(1)),
                       og.ogMap.values()))
        # ordering is important!
        assert map(itemgetter(0), og.ogMap['Anc1']) == ['HUMAN', 'CHIMP',
                                                        'BABOON']
        assert og.ogMap['Anc7'][0][0] == 'BABOON'
        assert og.ogMap['Anc7'][1][0] in ['CAT', 'DOG']
        assert og.ogMap['Anc7'][2][0] in ['CAT', 'DOG']

    def testMultipleOutgroupsOnRandomTrees(self):
        for tree in self.mcTrees:
            og = GreedyOutgroup()
            og.importTree(tree)
            og.greedy(candidateChildFrac=0.5, maxNumOutgroups=3)
            # make sure all entries have <= 3 outgroups.
            assert all(map(lambda x: len(x) <= 3, og.ogMap.values()))
            # and for all entries, the closest must be first.
            assert all(map(lambda x: x == sorted(x, key=itemgetter(1)),
                           og.ogMap.values()))

    def testDynamicOutgroupsOnRandomTrees(self):
        for tree, seqMap in zip(self.mcTrees, self.dummySeqMaps):
            degree = max([len(tree.getChildren(x)) for x in
                         tree.breadthFirstTraversal()])
            if degree < 8:
                og = DynamicOutgroup()
                og.edgeLen = 5
                og.importTree(tree, seqMap)
                og.compute(maxNumOutgroups=3)
                # make sure all entries have <= 3 outgroups.
                assert all(map(lambda x: len(x) <= 3, og.ogMap.values()))
                # and for all entries, the closest must be first.
                # (this will be true because all sequences are the same)
                assert all(map(lambda x: x == sorted(x, key=itemgetter(1)),
                               og.ogMap.values()))

    def testDynamicOutgroupsJustLeaves(self):
        tree = '((((HUMAN:0.006969,CHIMP:0.009727)Anc7:0.025291,BABOON:0.044568)Anc6:0.11,(MOUSE:0.072818,RAT:0.081244)Anc5:0.260342)Anc4:0.023260,((DOG:0.07,CAT:0.07)Anc3:0.087381,(PIG:0.06,COW:0.06)Anc2:0.104728)Anc1:0.04)Anc0;'
        mcTree = MultiCactusTree(NXNewick().parseString(tree, addImpliedRoots = False))
        mcTree.computeSubtreeRoots()
        og = DynamicOutgroup()
        og.importTree(mcTree, self.blanchetteSeqMap)
        og.compute(maxNumOutgroups=3, sequenceLossWeight=0.)
        # make sure all entries have <= 3 outgroups.
        assert all(map(lambda x: len(x) <= 3, og.ogMap.values()))
        # and for all entries, the closest must be first.
        assert all(map(lambda x: x == sorted(x, key=itemgetter(1)),
                       og.ogMap.values()))
        # ordering is important!
        assert og.ogMap['Anc1'][0][0] == 'HUMAN'
        assert og.ogMap['Anc7'][0][0] == 'BABOON'

        og = DynamicOutgroup()
        og.importTree(mcTree, self.blanchetteSeqMap)
        og.compute(maxNumOutgroups=3)
        # make sure all entries have <= 3 outgroups.
        assert all(map(lambda x: len(x) <= 3, og.ogMap.values()))

        # we keep dynamic outgroups sorted by distance too
        assert all(map(lambda x: x == sorted(x, key=itemgetter(1)),
                               og.ogMap.values()))
                        

    def testMultipleIdenticalRunsProduceSameResult(self):
        """The code now allows for multiple greedy() calls with different
        candidate sets, so that some outgroups can be 'preferred' over
        others without being the only candidates.
        Check that running greedy() multiple times with the same
        parameters gives the same result as running it once.
        """
        for tree in self.mcTrees:
            ogOnce = GreedyOutgroup()
            ogOnce.importTree(tree)
            ogOnce.greedy(maxNumOutgroups=3)
            ogMultipleTimes = GreedyOutgroup()
            ogMultipleTimes.importTree(tree)
            ogMultipleTimes.greedy(maxNumOutgroups=3)
            ogMultipleTimes.greedy(maxNumOutgroups=3)
            ogMultipleTimes.greedy(maxNumOutgroups=3)
            # make sure all entries have <= 3 outgroups.
            assert all(map(lambda x: len(x) <= 3, ogMultipleTimes.ogMap.values()))
            # and for all entries, the closest must be first.
            assert all(map(lambda x: x == sorted(x, key=itemgetter(1)),
                           ogMultipleTimes.ogMap.values()))
            # Check that the maps are equal. Can't compare them
            # directly since python will convert them to ordered
            # association lists.
            assert len(ogOnce.ogMap) == len(ogMultipleTimes.ogMap)
            for i in ogOnce.ogMap:
                assert i in ogMultipleTimes.ogMap
                assert ogOnce.ogMap[i] == ogMultipleTimes.ogMap[i]

    def testPreferredCandidateSets(self):
        """Test that running greedy() multiple times with different candidate
        sets will behave properly, i.e. keep all the existing outgroup
        assignments and fill in more on the second run."""
        for tree in self.mcTrees:
            ogOnce = GreedyOutgroup()
            ogOnce.importTree(tree)
            nodes = [j for j in tree.postOrderTraversal()]
            candidateSet = set([tree.getName(i) for i in random.sample(nodes, min(20, len(nodes)))])
            ogOnce.greedy(candidateSet=candidateSet, maxNumOutgroups=3)
            ogTwice = GreedyOutgroup()
            ogTwice.importTree(tree)
            ogTwice.greedy(candidateSet=candidateSet, maxNumOutgroups=3)
            ogTwice.greedy(maxNumOutgroups=3)
            # make sure all entries have <= 3 outgroups.
            assert all(map(lambda x: len(x) <= 3, ogTwice.ogMap.values()))
            # and for all entries, the closest must be first.
            assert all(map(lambda x: x == sorted(x, key=itemgetter(1)),
                           ogTwice.ogMap.values()))
            for node in ogTwice.ogMap:
                if node in ogOnce.ogMap:
                    # the ogMap entry in ogOnce should be a subset of the ogMap entry for ogTwice
                    oneRunOutgroups = ogOnce.ogMap[node]
                    twoRunOutgroups = ogTwice.ogMap[node]
                    assert len(twoRunOutgroups) >= len(oneRunOutgroups)
                    for i in oneRunOutgroups:
                        assert i in twoRunOutgroups

    def testNoOutgroupIsADescendantOfAnother(self):
        """No two outgroups should be on the same path to the root."""
        for tree in self.mcTrees:
            tree.nameUnlabeledInternalNodes()
            og = GreedyOutgroup()
            og.importTree(tree)
            og.greedy(maxNumOutgroups=3)
            for source in og.ogMap:
                for (sink1, _) in og.ogMap[source]:
                    for (sink2, _) in og.ogMap[source]:
                        if sink1 != sink2:
                            sink1Id = tree.nameToId[sink1]
                            sink2Id = tree.nameToId[sink2]
                            assert sink1Id not in tree.postOrderTraversal(sink2Id)
                            assert sink2Id not in tree.postOrderTraversal(sink1Id)
def main():
    unittest.main()

if __name__ == '__main__':
    main()
