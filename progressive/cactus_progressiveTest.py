#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Test the progressive workflow (completeley out of date right now)

"""
import unittest

import os
import sys
import random
import math
import copy

from optparse import OptionParser

from sonLib.bioio import getTempFile
from sonLib.bioio import getTempDirectory
from sonLib.bioio import TempFileTree
from sonLib.bioio import getRandomAlphaNumericString
from sonLib.bioio import printBinaryTree
from sonLib.bioio import newickTreeParser
from sonLib.bioio import system

from cladeDecomp import CladeDecomp
from cladeDecomp import getAllChildren
from cladeDecomp import subtree
from cladeDecomp import renameDbString

class TestCase(unittest.TestCase):
    def setUp(self):
        self.tempDir = getTempDirectory(os.getcwd())
        self.tempFiles = []
        self.fullTreeString = "(((A,B),(C,D)),((E,F),(G,H)))Root;"
        self.catTreeString = "(A,(B,(C,(D,(E,(F,(G,(H))))))));"
        self.sequences = ["A.fa", "B.fa", "C.fa", "D.fa", "E.fa", "F.fa", "G.fa", "H.fa"]
        self.optionsFull, args = OptionParser().parse_args()
        self.optionsFull.progWorkDir = self.tempDir
        self.optionsFull.speciesTree = self.fullTreeString
        self.optionsFull.cladeSize = 4
        self.optionsFull.cactusDiskDatabaseString = "<tokyo_cabinet database_dir=\"FOO\" />"
        self.optionsCat = copy.deepcopy(self.optionsFull)
        self.optionsCat.speciesTree = self.catTreeString
        self.optionsCat.cladeSize = 3
        unittest.TestCase.setUp(self)
        
        
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        system("rm -rf %s" % self.tempDir)
        for tempFile in self.tempFiles:
            os.remove(tempFile)
            
    def testNameInternalNodes(self):
        cdFull = CladeDecomp(self.optionsFull, self.sequences)
        namedFullTree = "(((A,B)Anc2,(C,D)Anc3)Anc0,((E,F)Anc4,(G,H)Anc5)Anc1)Root;"
        assert printBinaryTree(cdFull.tree, False) == namedFullTree
        cdCat = CladeDecomp(self.optionsCat, self.sequences)
        namedCatTree = "(A,(B,(C,(D,(E,(F,(G,H)Anc6)Anc5)Anc4)Anc3)Anc2)Anc1)Anc0;"
        assert printBinaryTree(cdCat.tree, False) == namedCatTree
    
    def testBuildSequenceIDLookup(self):

        cd = CladeDecomp(self.optionsFull, self.sequences)
        for i in ["A", "B", "C", "D", "E", "F", "G", "H"]:
            assert cd.lookup[i] == i + ".fa"
        for i in ["Anc0", "Anc1", "Anc2", "Anc3", "Anc4", "Anc5", "Root"]:
            assert cd.lookup[i] == cd.options.progWorkDir + "/" + i + ".fa"
    
    def testGetAllChildren(self):
        fullTree = newickTreeParser(self.fullTreeString)
        fullChildren = getAllChildren([fullTree.left, fullTree.right])
        assert fullTree.left.left in fullChildren
        assert fullTree.left.right in fullChildren
        assert fullTree.right.left in fullChildren
        assert fullTree.right.right in fullChildren
        assert len(fullChildren) is 4
        catTree = newickTreeParser(self.catTreeString)
        catChildren = getAllChildren([catTree])
        assert catTree.left in catChildren
        assert catTree.right in catChildren
        assert len(catChildren) is 2
    
    def testSubTree(self):
        cdFull = CladeDecomp(self.optionsFull, self.sequences)
        fullTree = cdFull.tree
        anc1 = fullTree.left
        anc2 = fullTree.right
        fullCladeTree = subtree(fullTree, [anc1, anc2])
        assert not fullCladeTree == fullTree
        assert fullCladeTree is not None
        assert printBinaryTree(fullCladeTree, False) == "(Anc0,Anc1)Root;"
        
        cdCat = CladeDecomp(self.optionsCat, self.sequences)
        catTree = cdCat.tree
        anc5 = catTree.right.right.right.right.right
        assert anc5.iD == "Anc5"
        catCladeTree = subtree(anc5, [anc5.left, anc5.right.left, anc5.right.right]) 
        assert not catCladeTree == anc5
        assert printBinaryTree(catCladeTree, False) == "(F,(G,H)Anc6)Anc5;"
    
    def testGetClade(self):
        cdFull = CladeDecomp(self.optionsFull, self.sequences)
        fullOptions, fullSequences, fullNodes = cdFull.getClade(cdFull.tree)
        ancList = ["Anc2", "Anc3", "Anc4", "Anc5"]
        ancList = map(lambda x: cdFull.options.progWorkDir + "/" + x +".fa", ancList)
        assert fullSequences == ancList
        treeNoDist = printBinaryTree(newickTreeParser(fullOptions.speciesTree), False)
        assert treeNoDist == "((Anc2,Anc3)Anc0,(Anc4,Anc5)Anc1)Root;" 
        assert fullOptions.cladeSize is self.optionsFull.cladeSize
        assert len(fullNodes) == cdFull.options.cladeSize
    
if __name__ == '__main__':
    unittest.main()
        
    