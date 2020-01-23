#!/usr/bin/env python3

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
"""
"""

import unittest
import os
import random
from operator import itemgetter
from sonLib.bioio import TestStatus
from sonLib.bioio import getTempDirectory
from sonLib.bioio import system

from cactus.progressive.multiCactusTree import MultiCactusTree
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
            if tree.size() < 50:
                mcTree = MultiCactusTree(tree)
                seqMap = dict()
                for i in mcTree.breadthFirstTraversal():
                    mcTree.setName(i, "Node%s" % str(i))
                    seqMap["Node%s" % str(i)] = self.tempFa
                mcTree.computeSubtreeRoots()
                mcTree.nameUnlabeledInternalNodes()
                self.mcTrees.append(mcTree)
                self.dummySeqMaps.append(seqMap)

        # Boreoeutherian tree
        borTree = '((((HUMAN:0.006969,CHIMP:0.009727)Anc7:0.025291,BABOON:0.044568)Anc6:0.11,(MOUSE:0.072818,RAT:0.081244)Anc5:0.260342)Anc4:0.023260,((DOG:0.07,CAT:0.07)Anc3:0.087381,(PIG:0.06,COW:0.06)Anc2:0.104728)Anc1:0.04)Anc0;'
        self.borMcTree = MultiCactusTree(NXNewick().parseString(borTree, addImpliedRoots=False))
        self.borMcTree.computeSubtreeRoots()
        self.borMcTree.nameUnlabeledInternalNodes()
        self.mcTrees.append(self.borMcTree)

        # Eutherian backbone tree
        backbone = '(((((((((((Homo_sapiens:0.00655,Pan_troglodytes:0.00684):0.00422,Gorilla_gorilla_gorilla:0.008964):0.009693,Pongo_abelii:0.01894):0.015511,Macaca_mulatta:0.043601):0.08444,Aotus_nancymaae:0.08):0.08,Microcebus_murinus:0.10612):0.043494,Galeopterus_variegatus:0.134937):0.04,((((Jaculus_jaculus:0.1,(Microtus_ochrogaster:0.14,(Mus_musculus:0.084509,Rattus_norvegicus:0.091589):0.047773):0.06015):0.122992,(Heterocephalus_glaber:0.1,(Cavia_porcellus:0.065629,(Chinchilla_lanigera:0.06,Octodon_degus:0.1):0.06):0.05):0.06015):0.05,Marmota_marmota:0.1):0.05,Oryctolagus_cuniculus:0.21569):0.04):0.040593,(((Sus_scrofa:0.12,(Orcinus_orca:0.069688,(Bos_taurus:0.04,Capra_hircus:0.04):0.09):0.045488):0.02,((Equus_caballus:0.109397,(Felis_catus:0.098612,(Canis_lupus_familiaris:0.052458,Mustela_putorius_furo:0.08):0.02):0.049845):0.02,(Pteropus_alecto:0.1,Eptesicus_fuscus:0.08):0.033706):0.03):0.025,Erinaceus_europaeus:0.278178):0.021227):0.023664,(((Loxodonta_africana:0.022242,Procavia_capensis:0.145358):0.076687,Chrysochloris_asiatica:0.04):0.05,Dasypus_novemcinctus:0.169809):0.02)backbone_root:0.234728,(Monodelphis_domestica:0.125686,Sarcophilus_harrisii:0.12):0.2151);'
        self.backboneTree = MultiCactusTree(NXNewick().parseString(backbone, addImpliedRoots=False))
        self.backboneTree.computeSubtreeRoots()
        self.backboneTree.nameUnlabeledInternalNodes()
        self.mcTrees.append(self.backboneTree)

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
        for event, seqLen in list(seqLens.items()):
            p = os.path.join(self.tempDir, event +".fa")
            with open(p, "w") as f:
                f.write(">%s\n" % event)
                f.write(''.join(['A'] * seqLen))
                f.write('\n')
            self.blanchetteSeqMap[event] = p

    def tearDown(self):
        unittest.TestCase.tearDown(self)
        system("rm -rf %s" % self.tempDir)

    @TestStatus.shortLength
    def testJustLeaves(self):
        og = GreedyOutgroup()
        og.importTree(self.borMcTree)
        candidates = set([self.borMcTree.getName(x) for x in self.borMcTree.getLeaves()])
        og.greedy(candidateSet=candidates, candidateChildFrac=2.)
        assert og.ogMap['Anc1'][0][0] == 'HUMAN'
        assert og.ogMap['Anc2'][0][0] in ['CAT', 'DOG']
        assert og.ogMap['Anc3'][0][0] in ['PIG', 'COW']
        assert og.ogMap['Anc4'][0][0] in ['CAT', 'DOG']
        assert og.ogMap['Anc5'][0][0] == 'HUMAN'
        assert og.ogMap['Anc6'][0][0] in ['CAT', 'DOG']
        assert og.ogMap['Anc7'][0][0] == 'BABOON'

    @TestStatus.shortLength
    def testHeightTable(self):
        """Make sure the height-table is calculated correctly."""
        og = GreedyOutgroup()
        og.importTree(self.borMcTree)
        htable = og.heightTable()
        self.assertEqual(htable[self.borMcTree.getNodeId('HUMAN')], 0)
        self.assertEqual(htable[self.borMcTree.getNodeId('PIG')], 0)
        self.assertEqual(htable[self.borMcTree.getNodeId('RAT')], 0)
        self.assertEqual(htable[self.borMcTree.getNodeId('Anc7')], 1)
        self.assertEqual(htable[self.borMcTree.getNodeId('Anc1')], 2)
        self.assertEqual(htable[self.borMcTree.getNodeId('Anc0')], 4)

    @TestStatus.shortLength
    def testZeroThreshold(self):
        """A threshold of 0 should produce outgroup sets that cause no additional depth in the resulting schedule."""
        tree = self.backboneTree
        og = GreedyOutgroup()
        og.importTree(tree)
        og.greedy(candidateSet=set(['Homo_sapiens', 'Mus_musculus']),threshold=0, maxNumOutgroups=3, candidateChildFrac=0.75)
        og.greedy(threshold=0, maxNumOutgroups=3, candidateChildFrac=0.75)
        htable = og.heightTable()
        for node, outgroups in list(og.ogMap.items()):
            for outgroup, _ in outgroups:
                # For the outgroup assignment to create no
                # additional dependencies, each outgroup must have
                # a height lower than the node it's outgroup to
                # (or be a leaf)
                self.assertTrue(htable[tree.getNodeId(outgroup)] < htable[tree.getNodeId(node)] \
                                or htable[tree.getNodeId(outgroup)] == 0)

    @TestStatus.shortLength
    def testCandidates(self):
        og = GreedyOutgroup()
        og.importTree(self.borMcTree)
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
        og.importTree(self.borMcTree)
        candidates = set(['HUMAN', 'CHIMP', 'RAT'])
        og.greedy(candidateSet=candidates, candidateChildFrac=1.0)
        assert og.ogMap['Anc1'][0][0] == 'Anc7'
        assert og.ogMap['Anc2'][0][0] == 'Anc7'
        assert og.ogMap['Anc3'][0][0] == 'Anc7'
        assert 'Anc4' not in og.ogMap
        assert og.ogMap['Anc5'][0][0] in ['HUMAN', 'CHIMP', 'Anc7']
        assert og.ogMap['Anc6'][0][0] == 'RAT'
        assert og.ogMap['Anc7'][0][0] == 'RAT'

    @TestStatus.shortLength
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

    @TestStatus.shortLength
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

    @TestStatus.shortLength
    def testMultipleOutgroups(self):
        og = GreedyOutgroup()
        og.importTree(self.borMcTree)
        og.greedy(candidateChildFrac=0.5, maxNumOutgroups=3)
        # make sure all entries have <= 3 outgroups.
        assert all([len(x) <= 3 for x in list(og.ogMap.values())])
        # and for all entries, the closest must be first.
        assert all([x == sorted(x, key=itemgetter(1)) for x in list(og.ogMap.values())])
        # ordering is important!
        assert list(map(itemgetter(0), og.ogMap['Anc4'])) == ['Anc1']
        assert list(map(itemgetter(0), og.ogMap['Anc7'])) == ['BABOON', 'Anc1',
                                                        'Anc5']
        # We avoid cycles, and choose post-order first, so this only
        # uses leaves.
        assert list(map(itemgetter(0), og.ogMap['Anc1'])) == ['HUMAN', 'CHIMP',
                                                        'BABOON']

    @TestStatus.shortLength
    def testMultipleOutgroupsJustLeaves(self):
        og = GreedyOutgroup()
        og.importTree(self.borMcTree)
        candidates = set([self.borMcTree.getName(x) for x in self.borMcTree.getLeaves()])
        og.greedy(candidateSet=candidates, candidateChildFrac=2.,
                  maxNumOutgroups=3)
        # make sure all entries have <= 3 outgroups.
        assert all([len(x) <= 3 for x in list(og.ogMap.values())])
        # and for all entries, the closest must be first.
        assert all([x == sorted(x, key=itemgetter(1)) for x in list(og.ogMap.values())])
        # ordering is important!
        assert list(map(itemgetter(0), og.ogMap['Anc1'])) == ['HUMAN', 'CHIMP',
                                                        'BABOON']
        assert og.ogMap['Anc7'][0][0] == 'BABOON'
        assert og.ogMap['Anc7'][1][0] in ['CAT', 'DOG']
        assert og.ogMap['Anc7'][2][0] in ['CAT', 'DOG']

    @TestStatus.shortLength
    def testMultipleOutgroupsOnRandomTrees(self):
        for tree in self.mcTrees:
            og = GreedyOutgroup()
            og.importTree(tree)
            og.greedy(candidateChildFrac=0.5, maxNumOutgroups=3)
            # make sure all entries have <= 3 outgroups.
            assert all([len(x) <= 3 for x in list(og.ogMap.values())])
            # and for all entries, the closest must be first.
            assert all([x == sorted(x, key=itemgetter(1)) for x in list(og.ogMap.values())])

    @TestStatus.shortLength
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
                assert all([len(x) <= 3 for x in list(og.ogMap.values())])
                # and for all entries, the closest must be first.
                # (this will be true because all sequences are the same)
                assert all([x == sorted(x, key=itemgetter(1)) for x in list(og.ogMap.values())])

    @TestStatus.shortLength
    def testDynamicOutgroupsJustLeaves(self):
        og = DynamicOutgroup()
        og.importTree(self.borMcTree, self.blanchetteSeqMap)
        og.compute(maxNumOutgroups=3, sequenceLossWeight=0.)
        # make sure all entries have <= 3 outgroups.
        assert all([len(x) <= 3 for x in list(og.ogMap.values())])
        # and for all entries, the closest must be first.
        assert all([x == sorted(x, key=itemgetter(1)) for x in list(og.ogMap.values())])
        # ordering is important!
        assert og.ogMap['Anc1'][0][0] == 'HUMAN'
        assert og.ogMap['Anc7'][0][0] == 'BABOON'

        og = DynamicOutgroup()
        og.importTree(self.borMcTree, self.blanchetteSeqMap)
        og.compute(maxNumOutgroups=3)
        # make sure all entries have <= 3 outgroups.
        assert all([len(x) <= 3 for x in list(og.ogMap.values())])

        # we keep dynamic outgroups sorted by distance too
        assert all([x == sorted(x, key=itemgetter(1)) for x in list(og.ogMap.values())])


    @TestStatus.shortLength
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
            assert all([len(x) <= 3 for x in list(ogMultipleTimes.ogMap.values())])
            # and for all entries, the closest must be first.
            assert all([x == sorted(x, key=itemgetter(1)) for x in list(ogMultipleTimes.ogMap.values())])
            # Check that the maps are equal. Can't compare them
            # directly since python will convert them to ordered
            # association lists.
            assert len(ogOnce.ogMap) == len(ogMultipleTimes.ogMap)
            for i in ogOnce.ogMap:
                assert i in ogMultipleTimes.ogMap
                assert ogOnce.ogMap[i] == ogMultipleTimes.ogMap[i]

    @TestStatus.shortLength
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
            assert all([len(x) <= 3 for x in list(ogTwice.ogMap.values())])
            # and for all entries, the closest must be first.
            assert all([x == sorted(x, key=itemgetter(1)) for x in list(ogTwice.ogMap.values())])
            for node in ogTwice.ogMap:
                if node in ogOnce.ogMap:
                    # the ogMap entry in ogOnce should be a subset of the ogMap entry for ogTwice
                    oneRunOutgroups = ogOnce.ogMap[node]
                    twoRunOutgroups = ogTwice.ogMap[node]
                    assert len(twoRunOutgroups) >= len(oneRunOutgroups)
                    for i in oneRunOutgroups:
                        assert i in twoRunOutgroups

    @TestStatus.shortLength
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
