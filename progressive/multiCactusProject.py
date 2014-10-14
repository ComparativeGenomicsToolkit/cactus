#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Basic interface to the multi cactus project xml file. 

"""
import unittest

import os
import xml.etree.ElementTree as ET
from xml.dom import minidom
import sys
import random
import math
import copy
import filecmp

from optparse import OptionParser

from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.experimentWrapper import ExperimentWrapper
from sonLib.nxnewick import NXNewick

class MultiCactusProject:
    def __init__(self):
        self.mcTree = None
        self.expMap = dict()
        self.inputSequences = []
        self.outputSequenceDir = None
        
    def readXML(self, path):
        xmlRoot = ET.parse(path).getroot()
        treeElem = xmlRoot.find("tree")
        self.mcTree = MultiCactusTree(NXNewick().parseString(treeElem.text, addImpliedRoots = False))
        self.expMap = dict()
        cactusPathElemList = xmlRoot.findall("cactus")
        for cactusPathElem in cactusPathElemList:
            nameElem = cactusPathElem.attrib["name"]
            pathElem = cactusPathElem.attrib["experiment_path"]
            self.expMap[nameElem] = pathElem
        self.inputSequences = xmlRoot.attrib["inputSequences"].split()
        self.outputSequenceDir = xmlRoot.attrib["outputSequenceDir"]
        self.mcTree.assignSubtreeRootNames(self.expMap)
        
    def writeXML(self, path):
        xmlRoot = ET.Element("multi_cactus")
        treeElem = ET.Element("tree")
        treeElem.text = NXNewick().writeString(self.mcTree)
        xmlRoot.append(treeElem)
        for name, expPath in self.expMap.items():
            cactusPathElem = ET.Element("cactus")
            cactusPathElem.attrib["name"] = name
            cactusPathElem.attrib["experiment_path"] = expPath
            xmlRoot.append(cactusPathElem)
        #We keep track of all the input sequences at the top level
        xmlRoot.attrib["inputSequences"] = " ".join(self.inputSequences)
        xmlRoot.attrib["outputSequenceDir"] = self.outputSequenceDir
        xmlFile = open(path, "w")
        xmlString = ET.tostring(xmlRoot)
        xmlString = minidom.parseString(xmlString).toprettyxml()
        xmlFile.write(xmlString)
        xmlFile.close()
    
    # find the sequence associated with an event name
    # by digging out the appropriate experiment file
    # doesn't work for the root!!!!
    def sequencePath(self, eventName):
        parentEvent = self.mcTree.getSubtreeRoot(eventName)           
        expPath = self.expMap[parentEvent]
        expElem = ET.parse(expPath).getroot()
        exp = ExperimentWrapper(expElem)
        seq = exp.getSequence(eventName)
        assert os.path.isfile(seq)
        return seq

    def getInputSequenceMap(self):
        """Return a map between event names and sequence paths.  Paths
        are different from above in that they are not taken from experiment
        xmls, but rather from directly from the project xml.
        """
        inputSequenceMap = dict()
        i = 0
        for node in self.mcTree.postOrderTraversal():
            if self.mcTree.isLeaf(node) is True:
                inputSequenceMap[self.mcTree.getName(node)] = \
                  self.inputSequences[i]
                i += 1
        assert i == len(self.inputSequences)
        return inputSequenceMap
        
    def getInputSequencePaths(self):
        """Get the set of input sequences for the multicactus tree
        """
        return self.inputSequences
    
    def getOutputSequenceDir(self):
        """The directory where the output sequences go
        """
        return self.outputSequenceDir
    
    def getConfigPath(self):
        return ExperimentWrapper(ET.parse(self.expMap.values()[0]).getroot()).getConfigPath()
      
if __name__ == '__main__':
    main()
        
    
