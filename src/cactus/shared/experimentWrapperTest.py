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
from sonLib.nxnewick import NXNewick

class TestCase(unittest.TestCase):
    
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.tree = '((((HUMAN:0.006969,CHIMP:0.009727):0.025291,BABOON:0.044568):0.11,(MOUSE:0.072818,RAT:0.081244):0.260342):0.02326,((DOG:0.07,CAT:0.07):0.087381,(PIG:0.06,COW:0.06):0.104728):0.04);'
        self.sequences = 'human.txt chimp.txt baboon.txt mouse.txt rat.txt dog.txt cat.txt pig.txt cow.txt' 
       
    def testOutgroups(self):
        xmlRoot = self.__makeXmlDummy(self.tree, self.sequences)
        exp = ExperimentWrapper(xmlRoot)
        assert NXNewick().writeString(exp.getTree()) == self.tree
        exp.addOutgroupSequence("outgroup", 1.3, "outgroup.fa")
        exp.addOutgroupSequence("outgroup2", 2.6, "outgroup2.fa")
        assert exp.getOutgroupEvents() == ["outgroup", "outgroup2"]
        seqMap = exp.buildSequenceMap()
        assert "outgroup" in seqMap
        assert seqMap["outgroup"] == "outgroup.fa"
        assert "outgroup2" in seqMap
        assert seqMap["outgroup2"] == "outgroup2.fa"
    def testSequenceMap(self):
        xmlRoot = self.__makeXmlDummy(self.tree, self.sequences)
        exp = ExperimentWrapper(xmlRoot)
        assert NXNewick().writeString(exp.getTree()) == self.tree
        
        seqMap = exp.buildSequenceMap()
        seqList = self.sequences.split()
        for i in seqList:
            assert seqMap[os.path.splitext(i)[0].upper()] == i
    
    def __makeXmlDummy(self, treeString, sequenceString):
        rootElem =  ET.Element("dummy")
        rootElem.attrib['species_tree'] = self.tree
        rootElem.attrib['sequences'] = self.sequences
        rootElem.append(self.__makeDiskElem())
        return rootElem
        
    def __makeDiskElem(self):
        diskElem = ET.Element("cactus_disk")
        confElem = ET.Element("st_kv_database_conf")
        confElem.attrib['type'] = 'redis'
        diskElem.append(confElem)
        dbElem = ET.Element('redis')
        confElem.append(dbElem)
        return diskElem 
            
    
def main():
    unittest.main()
    
if __name__ == '__main__':
    main()
