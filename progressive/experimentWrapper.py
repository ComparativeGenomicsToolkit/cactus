#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Interface to the cactus experiment xml file used
to read and motify an existing experiment

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
from sonLib.nxnewick import NXNewick
from cactus.shared.common import cactusRootPath

class ExperimentWrapper:
    def __init__(self, xmlRoot):
        self.xmlRoot = xmlRoot
        self.seqMap = self.buildSequenceMap()
        self.dbElem = self.getDbElem()

    def writeXML(self, path):
        xmlFile = open(path, "w")
        xmlString = ET.tostring(self.xmlRoot)
        xmlString = xmlString.replace("\n", "")
        xmlString = xmlString.replace("\t", "")
        xmlString = minidom.parseString(xmlString).toprettyxml()
        xmlFile.write(xmlString)
        xmlFile.close()

    def getDbDir(self):
        return self.dbElem.attrib["database_dir"]
    
    def setDbDir(self, path):
        self.dbElem.attrib["database_dir"] = path
        
    def getDbName(self):
        return self.dbElem.attrib["database_name"]
    
    def setDbName(self, name):
        if self.getDbType() == "kyoto_tycoon":
            assert os.path.splitext(name)[1] == ".kch"
        self.dbElem.attrib["database_name"] = name
    
    def getDbType(self):
        return self.dbElem.tag
    
    def getDbPort(self):
        assert self.getDbType() == "kyoto_tycoon"
        return int(self.dbElem.attrib["port"])
    
    def setDbPort(self, port):
        assert self.getDbType() == "kyoto_tycoon"
        self.dbElem.attrib["port"] = str(port)
    
    def getConfig(self):
        return self.xmlRoot.attrib["config"]
    
    def getTree(self):
        treeString = self.xmlRoot.attrib["species_tree"]
        return NXNewick().parseString(treeString)
    
    def getSequence(self, event):
        return self.seqMap[event]
    
    def getReferencePath(self):
        refElem = self.xmlRoot.find("reference")
        return refElem.attrib["path"]
    
    def setReferencePath(self, path):
        refElem = self.xmlRoot.find("reference")
        if refElem is None:
            refElem = ET.Element("reference")
            self.xmlRoot.append(refElem)
        refElem.attrib["path"] = path
    
    def getReferenceNameFromConfig(self):
        configElem = ET.parse(self.getConfig()).getroot()
        refElem = configElem.find("reference")
        return refElem.attrib["reference"]
        
    def getMAFPath(self):
        mafElem = self.xmlRoot.find("maf")
        return mafElem.attrib["path"]
    
    def setMAFPath(self, path):
        mafElem = self.xmlRoot.find("maf")
        if mafElem is None:
            mafElem = ET.Element("maf")
            self.xmlRoot.append(mafElem)
        mafElem.attrib["path"] = path
        
    def getOutgroupName(self):
        if self.xmlRoot.attrib.has_key("single_copy_species"):
            return self.xmlRoot.attrib["single_copy_species"]
        return None
    
    def addOutgroup(self, name, path, distance):
        assert self.getOutgroupName() is None
        tree = MultiCactusTree(self.getTree())
        tree.addOutgroup(name, distance)
        self.xmlRoot.attrib["species_tree"] = NXNewick().writeString(tree)   
        leaves = [tree.getName(i) for i in tree.getLeaves()]
        leaves.remove(name)     
        self.xmlRoot.attrib["required_species"] = ' '.join(leaves)
        self.xmlRoot.attrib["single_copy_species"] = name
        seqs = "%s %s"  % (self.xmlRoot.attrib["sequences"], path)
        self.xmlRoot.attrib["sequences"] = seqs
        self.seqMap[name] = path
    
    def getConfigPath(self):
        config = self.xmlRoot.attrib["config"]
        if config == 'default':
            dir = os.path.join(cactusRootPath(), "pipeline")
            config = os.path.join(dir, "cactus_workflow_config.xml")
        if config == 'defaultProgressive':
            dir = os.path.join(cactusRootPath(), "progressive")
            config = os.path.join(dir, "cactus_progressive_workflow_config.xml")
        return config
    
    def setConfigPath(self, path):
        self.xmlRoot.attrib["config"] = path
        
    def getDiskDatabaseString(self):
        conf = self.xmlRoot.find("cactus_disk").find("st_kv_database_conf")
        return ET.tostring(conf).replace("\n", "").replace("\t","")
    
    # map event names to sequence paths
    def buildSequenceMap(self):
        tree = self.getTree()
        sequenceString = self.xmlRoot.attrib["sequences"]
        sequences = sequenceString.split()
        nameIterator = iter(sequences)
        seqMap = dict()
        for node in tree.postOrderTraversal():
            if tree.isLeaf(node):
                seqMap[tree.getName(node)] = nameIterator.next()
        return seqMap
    
    def getDbElem(self):
        diskElem = self.xmlRoot.find("cactus_disk")
        confElem = diskElem.find("st_kv_database_conf")
        typeString = confElem.attrib["type"]
        return confElem.find(typeString)  

    # load in a new tree (using input seqMap if specified,
    # current one otherwise
    def updateTree(self, tree, seqMap = None):
        if seqMap is not None:
            self.seqMap = seqMap
        newMap = dict()
        treeString = NXNewick().writeString(tree)
        self.xmlRoot.attrib["species_tree"] = treeString
        sequences = "" 
        for node in tree.postOrderTraversal():
            if tree.isLeaf(node):
                nodeName = tree.getName(node)
                if len(sequences) > 0:
                    sequences += " "                  
                sequences += seqMap[nodeName]
                newMap[nodeName] = seqMap[nodeName]
        self.xmlRoot.attrib["sequences"] = sequences
        self.seqMap = newMap
        
        
    