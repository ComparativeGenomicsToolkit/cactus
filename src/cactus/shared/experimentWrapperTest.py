#!/usr/bin/env python3

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
"""
"""

import unittest
import os
import xml.etree.ElementTree as ET
from sonLib.bioio import TestStatus
from cactus.shared.experimentWrapper import ExperimentWrapper
from sonLib.nxnewick import NXNewick

class TestCase(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)
        self.tree = NXNewick().parseString('((((HUMAN:0.006969,CHIMP:0.009727)anc2:0.025291,BABOON:0.044568)anc1:0.11,(MOUSE:0.072818,RAT:0.081244):0.260342):0.02326,((DOG:0.07,CAT:0.07):0.087381,(PIG:0.06,COW:0.06):0.104728):0.04);')
        self.xmlRoot = self.__makeXmlDummy()
        self.exp = ExperimentWrapper(self.xmlRoot)
        self.exp.setTree(self.tree)
        self.seqMap = {'HUMAN': 'human.txt',
                       'CHIMP': 'chimp.txt',
                       'BABOON': 'baboon.txt',
                       'MOUSE': 'mouse.txt',
                       'RAT': 'rat.txt',
                       'DOG': 'dog.txt',
                       'CAT': 'cat.txt',
                       'PIG': 'pig.txt',
                       'COW': 'cow.txt'}
        self.exp.setRootGenome('anc1')
        self.exp.setRootReconstructed(True)
        self.exp.setOutgroupGenomes(['MOUSE', 'RAT', 'DOG', 'CAT', 'PIG', 'COW'])
        for genome, seq in list(self.seqMap.items()):
            # These aren't real IDs, but should still work for our
            # purposes
            self.exp.setSequenceID(genome, seq)

    @TestStatus.shortLength
    def testGetSequencePath(self):
        for genome, seq in list(self.seqMap.items()):
            self.assertEqual(self.exp.getSequenceID(genome), seq)

        # Should not be any entries for genomes not in the tree
        self.assertEqual(self.exp.getSequenceID('DUCK'), None)

    @TestStatus.shortLength
    def testChangingSequencePaths(self):
        """Tests that changing a sequence path persists correctly."""
        self.exp.setSequenceID('HUMAN', 'human2.txt')
        self.assertEqual(self.exp.getSequenceID('HUMAN'), 'human2.txt')
        # Reload the wrapper and try again
        self.exp = ExperimentWrapper(self.xmlRoot)
        self.assertEqual(self.exp.getSequenceID('HUMAN'), 'human2.txt')

    @TestStatus.shortLength
    def testOutgroups(self):
        self.assertEqual(set(self.exp.getOutgroupGenomes()),
                         set(['MOUSE', 'RAT', 'DOG', 'CAT', 'PIG', 'COW']))
        self.exp.setOutgroupGenomes([])
        self.assertEqual(self.exp.getOutgroupGenomes(), [])

    @TestStatus.shortLength
    def testRootGenome(self):
        self.assertEqual(self.exp.getRootGenome(), 'anc1')
        self.exp.setRootGenome('anc2')
        self.assertEqual(self.exp.getRootGenome(), 'anc2')

    @TestStatus.shortLength
    def testSetTree(self):
        # A modfied version, with fewer genomes and a new one
        tree2 = NXNewick().parseString('((HUMAN:0.006969,CHIMP:0.009727):0.025291,BABOON:0.044568,ARMADILLO:1.0);')
        self.exp.setTree(tree2)
        self.assertEqual(set(self.exp.getGenomesWithSequence()), set(['HUMAN', 'CHIMP', 'BABOON']))

    def __makeXmlDummy(self):

        rootElem =  ET.Element("dummy")
        rootElem.append(self.__makeDiskElem())
        return rootElem

    def __makeDiskElem(self):
        diskElem = ET.Element("cactus_disk")
        confElem = ET.Element("st_kv_database_conf")
        confElem.attrib['type'] = 'kyoto_tycoon'
        diskElem.append(confElem)
        dbElem = ET.Element('kyoto_tycoon')
        confElem.append(dbElem)
        return diskElem


def main():
    unittest.main()

if __name__ == '__main__':
    main()
