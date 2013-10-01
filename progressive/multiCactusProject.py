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
        
    def readXML(self, path):
        xmlRoot = ET.parse(path).getroot()
        treeElem = xmlRoot.find("tree")
        self.mcTree = MultiCactusTree(NXNewick().parseString(treeElem.text))
        self.expMap = dict()
        cactusPathElemList = xmlRoot.findall("cactus")
        for cactusPathElem in cactusPathElemList:
            nameElem = cactusPathElem.attrib["name"]
            pathElem = cactusPathElem.attrib["experiment_path"]
            self.expMap[nameElem] = pathElem
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
        seq = exp.getOutputSequence(eventName)
        assert os.path.isfile(seq)
        return seq
    
    def getInputSequencePaths(self):
        """Get the set of input sequences for the multicactus tree
        """
        sequences = set()
        for expPath in self.expMap.values():
            eW = ExperimentWrapper(ET.parse(expPath).getroot())
            for i in eW.getSequences():
                sequences.add(i)
        return list(sequences)

    def getAnExperimentWrapper(self):
        #Load an experiment file to get stuff that is invariant between experiments
        return ExperimentWrapper(ET.parse(self.expMap.values()[0]).getroot())
    
    def getConfigPath(self):
        return self.getAnExperimentWrapper().getConfigPath()
      
if __name__ == '__main__':
    main()
        
    