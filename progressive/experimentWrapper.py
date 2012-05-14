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
from cactus.progressive.configWrapper import ConfigWrapper
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
        if "database_dir" in self.dbElem.attrib:
            dbDir = self.dbElem.attrib["database_dir"]
            if len(dbDir) > 0:
                return self.dbElem.attrib["database_dir"]
        return None
    
    def setDbDir(self, path):
        self.dbElem.attrib["database_dir"] = path
        
    def getDbName(self):
        if self.getDbType() == "kyoto_tycoon" and self.getDbInMemory() == True:
            return ""
        return self.dbElem.attrib["database_name"]
    
    def setDbName(self, name):
        if self.getDbType() == "kyoto_tycoon":
            if self.getDbInMemory() == True:
                if "database_name" in self.dbElem.attrib:
                    del self.dbElem.attrib["database_name"]
                return
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
    
    def getDbHost(self):
        assert self.getDbType() == "kyoto_tycoon"
        if "host" in self.dbElem.attrib:
            return self.dbElem.attrib["host"]
        return None
        
    def getDbServerOptions(self):
        assert self.getDbType() == "kyoto_tycoon"
        if "server_options" in self.dbElem.attrib:
            return self.dbElem.attrib["server_options"]
        return None
    
    def getDbTuningOptions(self):
        assert self.getDbType() == "kyoto_tycoon"
        if "tuning_options" in self.dbElem.attrib:
            return self.dbElem.attrib["tuning_options"]
        return None

    def getDbCreateTuningOptions(self):
        assert self.getDbType() == "kyoto_tycoon"
        if "create_tuning_options" in self.dbElem.attrib:
            return self.dbElem.attrib["create_tuning_options"]
        return None

    def getDbReadTuningOptions(self):
        assert self.getDbType() == "kyoto_tycoon"
        if "read_tuning_options" in self.dbElem.attrib:
            return self.dbElem.attrib["read_tuning_options"]
        return None
    
    def getDbInMemory(self):
        assert self.getDbType() == "kyoto_tycoon"
        if "in_memory" in self.dbElem.attrib:
            val = self.dbElem.attrib["in_memory"]
            retVal = val.lower() == "true" or val == "1"
            assert (not retVal or "database_name" not in self.dbElem.attrib)
            return retVal
        return False
    
    def getDbSnapshot(self):
        assert self.getDbType() == "kyoto_tycoon"
        if "snapshot" in self.dbElem.attrib:
            val = self.dbElem.attrib["snapshot"]
            return val.lower() == "true" or val == "1"
        return self.getDbInMemory()
    
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
        
    def getHALPath(self):
        halElem = self.xmlRoot.find("hal")
        if halElem is None:
            return None
        return halElem.attrib["path"]
    
    def setHALPath(self, path):
        halElem = self.xmlRoot.find("hal")
        if halElem is None:
            halElem = ET.Element("hal")
            self.xmlRoot.append(halElem)
        assert os.path.splitext(path)[1] == ".hal"
        halElem.attrib["path"] = path
        
    def getMAFPath(self):
        halElem = self.xmlRoot.find("hal")
        if halElem is None:
            return None
        return halElem.attrib["path"]
    
    def setMAFPath(self, path):
        halElem = self.xmlRoot.find("hal")
        if halElem is None:
            halElem = ET.Element("hal")
            self.xmlRoot.append(halElem)
        assert os.path.splitext(path)[1] == ".maf"
        halElem.attrib["path"] = path
        
    def getJoinMAF(self):
        configElem = ET.parse(self.getConfig()).getroot()
        config = ConfigWrapper(configElem)
        return config.getJoinMaf()
            
    def getOutgroupEvents(self):
        if self.xmlRoot.attrib.has_key("outgroup_events"):
            return self.xmlRoot.attrib["outgroup_events"].split()
        return []
        
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
    
    # sets the outgroup.
    def setOutgroup(self, ogName, ogDist, ogPath):
        tree = MultiCactusTree(self.getTree())
        leaves = [tree.getName(i) for i in tree.getLeaves()]
        
        if ogName is not None:
            tree.addOutgroup(ogName, ogDist)
            self.xmlRoot.attrib["species_tree"] = NXNewick().writeString(tree)   
            seqs = "%s %s"  % (self.xmlRoot.attrib["sequences"], ogPath)
            self.xmlRoot.attrib["sequences"] = seqs
            self.seqMap[ogName] = ogPath
            self.xmlRoot.attrib["outgroup_events"] = ogName
      