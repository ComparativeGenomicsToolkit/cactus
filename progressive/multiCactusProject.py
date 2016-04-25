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

from toil.lib.bioio import logger

class MultiCactusProject:
    def __init__(self):
        self.mcTree = None
        self.expMap = dict()
        self.expIDMap = None
        self.inputSequences = []
        self.outputSequenceDir = None
        self.inputSequenceIDs = None
        self.outputSequenceIDMap = None
        self.configID = None
        
    def readXML(self, path, jobStore=None):
        xmlRoot = ET.parse(path).getroot()
        treeElem = xmlRoot.find("tree")
        self.mcTree = MultiCactusTree(NXNewick().parseString(treeElem.text, addImpliedRoots = False))
        self.expMap = dict()
        self.expIDMap = dict()
        cactusPathElemList = xmlRoot.findall("cactus")
        for cactusPathElem in cactusPathElemList:
            nameElem = cactusPathElem.attrib["name"]
            pathElem = cactusPathElem.attrib["experiment_path"]
            self.expMap[nameElem] = pathElem
            if "experiment_id" in cactusPathElem.attrib:
                self.expIDMap[nameElem] =cactusPathElem.attrib["experiment_id"]
            elif jobStore:
                self.expIDMap[nameElem] = jobStore.importFile("file://" + pathElem)
        self.inputSequences = xmlRoot.attrib["inputSequences"].split()
        if "inputSequenceIDs" in xmlRoot.attrib:
            self.inputSequenceIDs = xmlRoot.attrib["inputSequenceIDs"].split()
        if "outputSequenceIDs" in xmlRoot.attrib:
            self.outputSequenceIDMap = dict(zip(xmlRoot.attrib["outputSequenceIDs"].split(),
                                                xmlRoot.attrib["outputSequenceNames"].split()))

        logger.info("xmlRoot = %s" % ET.tostring(xmlRoot))
        if "configID" in xmlRoot.attrib:
            self.configID = xmlRoot.attrib["configID"]
            
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
            if self.expIDMap:
                cactusPathElem.attrib["experiment_id"] = self.expIDMap[name]
            xmlRoot.append(cactusPathElem)
        #We keep track of all the input sequences at the top level
        xmlRoot.attrib["inputSequences"] = " ".join(self.inputSequences)
        xmlRoot.attrib["outputSequenceDir"] = self.outputSequenceDir
        if self.inputSequenceIDs:
            xmlRoot.attrib["inputSequenceIDs"] = " ".join(self.inputSequenceIDs)
        if self.outputSequenceIDMap:
            xmlRoot.attrib["outputSequenceIDs"] = " ".join(self.outputSequenceIDMap.values())
            xmlRoot.attrib["outputSequenceNames"] = " ".join(self.outputSequenceIDMap.keys())
        if self.configID:
            xmlRoot.attrib["configID"] = self.configID

        xmlFile = open(path, "w")
        xmlString = ET.tostring(xmlRoot)
        xmlString = minidom.parseString(xmlString).toprettyxml()
        xmlFile.write(xmlString)
        xmlFile.close()

    def writeXMLToFileStore(self, fileStore):
        tmpXML = os.path.join(fileStore.getLocalTempDir(), "tmpXML")
        self.writeXML(tmpXML)
        return fileStore.writeGlobalFile(tmpXML)

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



    def getInputSequenceIDs(self):
        """Get the set of input sequences for the multicactus tree
        """
        return self.inputSequenceIDs
    
    def getInputSequencePaths(self):
        return self.inputSequences
    
    def getOutputSequenceDir(self):
        """The directory where the output sequences go
        """
        return self.outputSequenceDir
        
    def setOutputSequenceIDs(self, outputSequenceIDs):
        self.outputSequenceIDMap = dict()
        i = 0
        for node in self.mcTree.postOrderTraversal():
            if self.mcTree.isLeaf(node) is True:
                self.outputSequenceIDMap[self.mcTree.getName(node)] = \
                  outputSequenceIDs[i]
                i += 1
        assert i == len(outputSequenceIDs)
        
    def getOutputSequenceIDs(self):
        return self.outputSequenceIDs
    
    def getConfigPath(self):
        return ExperimentWrapper(ET.parse(self.expMap.values()[0]).getroot()).getConfigPath()

    def setConfigID(self, configID):
        self.configID = configID

    def getConfigID(self):
        return self.configID

    def setInputSequenceIDs(self, inputSequenceIDs):
        self.inputSequenceIDs = inputSequenceIDs
      
if __name__ == '__main__':
    main()
        
    
