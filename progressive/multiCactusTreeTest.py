#!/usr/bin/env python

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
    
    def testSanity(self):
        parser = NXNewick()
        mcTree1 = MultiCactusTree(parser.parseString(self.tree1))
        tree1String = NXNewick().writeString(mcTree1)  
        assert tree1String == self.tree1
        mcTree2 = MultiCactusTree(parser.parseString(self.tree2), subtreeSize = 3)
        tree2String = NXNewick().writeString(mcTree2)
        assert tree2String == self.tree2
        
    def testSubtrees(self):
        roots1 = ["Anc0", "Anc1", "Anc2", "Anc3", "Anc4", "Anc5", "Anc6", "Anc7"]
        roots1a = ["Anc0", "Anc2", "Anc3", "Anc5", "Anc6", "Anc7"]
        roots2 = ["Anc0", "Anc1", "Anc2", "Anc3", "Anc4", "Anc5"]
        
        subTree1_a6 = '(Anc7:0.025291,BABOON:0.044568)Anc6;'
        subTree1a_a0 = '((Anc6:0.11,Anc5:0.260342)Anc4:0.02326,(Anc3:0.087381,Anc2:0.104728)Anc1:0.04)Anc0;'
        subTree2_a3 = '(monkey:100.8593,cat:47.14069)Anc3;'
        
        trueRoots = [roots1, roots1a, roots2]
        trueSubtrees = [subTree1_a6, subTree1a_a0, subTree2_a3]
        trees = [self.mcTree1, self.mcTree1a, self.mcTree2]
        ancs = ["Anc6", "Anc0", "Anc3"]
        
        for i in range(0, 3):
            roots = trees[i].getSubtreeRootNames()
            assert sorted(roots) == sorted(trueRoots[i])
            subtree = trees[i].extractSubTree(ancs[i])
            subtree = NXNewick().writeString(subtree)
            assert subtree == trueSubtrees[i]        
    
    def testAddSelf(self):
        trueSelf = '((((((COW:0.06)COW_self:0.06,(PIG:0.06)PIG_self:0.06)Anc2:0.104728)Anc2_self:0.104728,(((CAT:0.07)CAT_self:0.07,(DOG:0.07)DOG_self:0.07)Anc3:0.087381)Anc3_self:0.087381)Anc1)Anc1_self:0.04,(((((RAT:0.081244)RAT_self:0.081244,(MOUSE:0.072818)MOUSE_self:0.072818)Anc5:0.260342)Anc5_self:0.260342,(((BABOON:0.044568)BABOON_self:0.044568,(((CHIMP:0.009727)CHIMP_self:0.009727,(HUMAN:0.006969)HUMAN_self:0.006969)Anc7:0.025291)Anc7_self:0.025291)Anc6:0.11)Anc6_self:0.11)Anc4)Anc4_self:0.02326)Anc0;'
        tree = MultiCactusTree(self.mcTree1)
        tree.nameUnlabeledInternalNodes()
        tree.computeSubtreeRoots()
        tree.addSelfEdges()
        treeString = NXNewick().writeString(tree)
        assert treeString == trueSelf
    
    def testAddOutgroup(self):
        trueOg = '((((HUMAN:0.006969,CHIMP:0.009727)Anc7:0.025291,BABOON:0.044568)Anc6:0.11,(MOUSE:0.072818,RAT:0.081244)Anc5:0.260342)Anc4:0.02326,((DOG:0.07,CAT:0.07)Anc3:0.087381,(PIG:0.06,COW:0.06)Anc2:0.104728)Anc1:0.04,outgroup:1.7)Anc0;'
        tree = MultiCactusTree(self.mcTree1)
        tree.nameUnlabeledInternalNodes()
        tree.computeSubtreeRoots()
        tree.addOutgroup("outgroup", 1.7)        
        treeString = NXNewick().writeString(tree)
        assert treeString == trueOg
        
        trueLeafOg = "(A:1.1,outgroup:1.1);"
        leafTreeString = "A;"
        parser = NXNewick()
        leafTree = MultiCactusTree(parser.parseString(leafTreeString))
        leafTree.nameUnlabeledInternalNodes()
        leafTree.computeSubtreeRoots()
        leafTree.addOutgroup("outgroup", 2.2)
        leafTreeOutString = NXNewick().writeString(leafTree)
        assert leafTreeOutString == trueLeafOg

    def testMakeSelfName(self):
        tree = MultiCactusTree()
        name = "chimp.chr6"
        selfname = "chimp_self.chr6"
        outname = tree.makeSelfName(name, "_self")
        assert outname == selfname
                   
    def __generateTrees(self):
        self.tree1 = '((((HUMAN:0.006969,CHIMP:0.009727):0.025291,BABOON:0.044568):0.11,(MOUSE:0.072818,RAT:0.081244):0.260342):0.02326,((DOG:0.07,CAT:0.07):0.087381,(PIG:0.06,COW:0.06):0.104728):0.04);'
        self.tree2 = '((raccoon:19.19959,bear:6.80041):0.846,((sea_lion:11.997,seal:12.003):7.52973,((monkey:100.8593,cat:47.14069):20.59201,weasel:18.87953):2.0946):3.87382,dog:25.46154);'
        parser = NXNewick()
        self.mcTree1 = MultiCactusTree(parser.parseString(self.tree1))
        self.mcTree1a = MultiCactusTree(parser.parseString(self.tree1), subtreeSize = 4)
        self.mcTree2 = MultiCactusTree(parser.parseString(self.tree2), subtreeSize = 3)
        self.mcTree1.nameUnlabeledInternalNodes()
        self.mcTree1a.nameUnlabeledInternalNodes()
        self.mcTree2.nameUnlabeledInternalNodes()
        self.mcTree1.computeSubtreeRoots()
        self.mcTree1a.computeSubtreeRoots()
        self.mcTree2.computeSubtreeRoots()
        
    
def main():
    unittest.main()
    
if __name__ == '__main__':
    main()
