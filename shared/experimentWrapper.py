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
from sonLib.bioio import getRandomAlphaNumericString
from sonLib.bioio import system

from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.configWrapper import ConfigWrapper
from sonLib.nxnewick import NXNewick
from cactus.shared.common import cactusRootPath

class DbElemWrapper(object):
    def __init__(self, confElem):
        typeString = confElem.attrib["type"]
        dbElem = confElem.find(typeString)
        self.dbElem = dbElem
        self.confElem = confElem
        
    def check(self):
        """Function checks the database conf is as expected and creates useful exceptions
        if not"""
        dataString = ET.tostring(self.confElem)
        if self.confElem.tag != "st_kv_database_conf":
            raise RuntimeError("The database conf string is improperly formatted: %s" % dataString)
        if not self.confElem.attrib.has_key("type"):
            raise RuntimeError("The database conf string does not have a type attrib: %s" % dataString)
        typeString = self.confElem.attrib["type"]
        if typeString == "tokyo_cabinet":
            tokyoCabinet = self.confElem.find("tokyo_cabinet")
            if tokyoCabinet == None:
                raise RuntimeError("Database conf is of type tokyo cabinet but there is no nested tokyo cabinet tag: %s" % dataString)
            if not tokyoCabinet.attrib.has_key("database_dir"):
                raise RuntimeError("The tokyo cabinet tag has no database_dir tag: %s" % dataString)
        elif typeString == "kyoto_tycoon":
            kyotoTycoon = self.confElem.find("kyoto_tycoon")
            if kyotoTycoon == None:
                raise RuntimeError("Database conf is of kyoto tycoon but there is no nested kyoto tycoon tag: %s" % dataString)
            if not set(("host", "port", "database_dir")).issubset(set(kyotoTycoon.attrib.keys())):
                raise RuntimeError("The kyoto tycoon tag has a missing attribute: %s" % dataString)
        else:
            raise RuntimeError("Unrecognised database type in conf string: %s" % typeString)

    def getDbElem(self):
        return self.dbElem

    def getConfString(self):
        return ET.tostring(self.confElem)

    def getDbDir(self): #Replacement for getDatabaseString
        if "database_dir" in self.dbElem.attrib:
            dbDir = self.dbElem.attrib["database_dir"]
            if len(dbDir) > 0:
                if dbDir[-1] == '/' and lendbDir > 1:
                    dbDir = dbDir[:-1]
                return dbDir
        return None
    
    def setDbDir(self, path):
        if path[-1] == '/' and len(path) > 1:
            self.dbElem.attrib["database_dir"] = path[:-1]
        else:
            self.dbElem.attrib["database_dir"] = path
        
    def getDbName(self):
        #if self.getDbType() == "kyoto_tycoon" and self.getDbInMemory() == True:
        #    return ""
        if "database_name" not in self.dbElem.attrib:
            return None
        return self.dbElem.attrib["database_name"]
    
    def setDbName(self, name):
        #if self.getDbType() == "kyoto_tycoon":
        #    if self.getDbInMemory() == True:
        #        if "database_name" in self.dbElem.attrib:
        #            del self.dbElem.attrib["database_name"]
        #        return
        #    assert os.path.splitext(name)[1] == ".kch"            
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
    
    def setDbHost(self, host):
        assert self.getDbType() == "kyoto_tycoon"
        self.dbElem.attrib["host"] = host
        
    def getDbServerOptions(self):
        assert self.getDbType() == "kyoto_tycoon"
        if "server_options" in self.dbElem.attrib:
            return self.dbElem.attrib["server_options"]
        return None

    def setDbServerOptions(self, options):
        assert self.getDbType() == "kyoto_tycoon"
        self.dbElem.attrib["server_options"] = str(options)
    
    def getDbTuningOptions(self):
        assert self.getDbType() == "kyoto_tycoon"
        if "tuning_options" in self.dbElem.attrib:
            return self.dbElem.attrib["tuning_options"]
        return None

    def setDbTuningOptions(self, options):
        assert self.getDbType() == "kyoto_tycoon"
        self.dbElem.attrib["tuning_options"] = str(options)

    def getDbCreateTuningOptions(self):
        assert self.getDbType() == "kyoto_tycoon"
        if "create_tuning_options" in self.dbElem.attrib:
            return self.dbElem.attrib["create_tuning_options"]
        return None

    def setDbCreateTuningOptions(self, options):
        assert self.getDbType() == "kyoto_tycoon"
        self.dbElem.attrib["create_tuning_options"] = str(options)

    def getDbReadTuningOptions(self):
        assert self.getDbType() == "kyoto_tycoon"
        if "read_tuning_options" in self.dbElem.attrib:
            return self.dbElem.attrib["read_tuning_options"]
        return None

    def setDbReadTuningOptions(self, options):
        assert self.getDbType() == "kyoto_tycoon"
        self.dbElem.attrib["read_tuning_options"] = str(options)

    def getDbInMemory(self):
        assert self.getDbType() == "kyoto_tycoon"
        if "in_memory" in self.dbElem.attrib:
            val = self.dbElem.attrib["in_memory"]
            retVal = val.lower() == "true" or val == "1"
            #assert (not retVal or "database_name" not in self.dbElem.attrib)
            return retVal
        return False

    def setDbInMemory(self, inMemory):
        assert self.getDbType() == "kyoto_tycoon"
        self.dbElem.attrib["in_memory"] = str(int(inMemory))
    
    def getDbSnapshot(self):
        assert self.getDbType() == "kyoto_tycoon"
        if "snapshot" in self.dbElem.attrib:
            val = self.dbElem.attrib["snapshot"]
            return val.lower() == "true" or val == "1"
        return self.getDbInMemory()

    def setDbSnapshot(self, snapshot):
        assert self.getDbType() == "kyoto_tycoon"
        self.dbElem.attrib["snapshot"] = str(int(snapshot))
    
    def cleanupDb(self): #Replacement for cleanupDatabase
        """Removes the database that was created.
        """
        if self.getDbType() == "kyoto_tycoon":
            system("ktremotemgr clear -port %s -host %s" % (self.getDbPort(), self.getDbHost()))
            system("rm -rf %s" % self.getDbDir())
        else:
            assert self.getDbDir() != None
            system("rm -rf %s" % self.getDbDir())

class ExperimentWrapper(DbElemWrapper):
    def __init__(self, xmlRoot):
        self.diskElem = xmlRoot.find("cactus_disk")
        confElem = self.diskElem.find("st_kv_database_conf")
        super(ExperimentWrapper, self).__init__(confElem)
        self.xmlRoot = xmlRoot
        self.seqMap = self.buildSequenceMap()
        
    @staticmethod
    def createExperimentWrapper(sequences, newickTreeString, outputDir,
                 outgroupEvents=None,
                 databaseConf=None, configFile=None, 
                 halFile=None, fastaFile=None,
                 constraints=None, progressive=False,
                 outputSequenceDir=None):
        #Establish the basics
        rootElem =  ET.Element("cactus_workflow_experiment")
        rootElem.attrib['species_tree'] = newickTreeString
        rootElem.attrib['sequences'] = " ".join(sequences)
        #Stuff for the database
        database = ET.SubElement(rootElem, "cactus_disk")
        if databaseConf != None:
            database.append(databaseConf)
        else:
            databaseConf = ET.SubElement(database, "st_kv_database_conf")
            databaseConf.attrib["type"] = "tokyo_cabinet"
            tokyoCabinet = ET.SubElement(databaseConf, "tokyo_cabinet")
        self = ExperimentWrapper(rootElem) 
        #Output dir
        self.setOutputDir(outputDir)
        #Database name
        if self.getDbName() == None:
            self.setDbName("cactusDisk_%s_%i" % (getRandomAlphaNumericString(), os.getpid())) #Needs to be unique
        #Database dir
        if self.getDbDir() == None:
            self.setDbDir(os.path.join(outputDir, self.getDbName()))
        #Hal file
        if halFile != None:
            self.setHALPath(halFile) 
        #Fasta file
        if fastaFile != None:
            self.setHALFastaPath(fastaFile)
        #Setup the config
        self.setConfigPath("default")
        if progressive == True:
            self.setConfigPath("defaultProgressive")
        if configFile != None:
            self.setConfigPath(configFile)
        #Constraints
        if constraints != None:
            self.setConstraintsFilePath(constraints)
        #Outgroup events
        if outgroupEvents != None:
            self.setOutgroupEvents(outgroupEvents)
        #Output sequences
        if outputSequenceDir != None:
            self.setOutputSequenceDir(outputSequenceDir)
        return self

    def writeXML(self, path): #Replacement for writeExperimentFile
        xmlFile = open(path, "w")
        xmlString = ET.tostring(self.xmlRoot)
        xmlString = xmlString.replace("\n", "")
        xmlString = xmlString.replace("\t", "")
        xmlString = minidom.parseString(xmlString).toprettyxml()
        xmlFile.write(xmlString)
        xmlFile.close()
    
    def getConfig(self):
        return self.xmlRoot.attrib["config"]
    
    def getTree(self):
        treeString = self.xmlRoot.attrib["species_tree"]
        return NXNewick().parseString(treeString, addImpliedRoots = False)
    
    def setSequences(self, sequences):
        self.xmlRoot.attrib["sequences"] = " ".join(sequences)
        self.seqMap = self.buildSequenceMap()
    
    def getSequences(self):
        return self.xmlRoot.attrib["sequences"].split()
    
    def getSequence(self, event):
        return self.seqMap[event]
    
    def getReferencePath(self):
        refElem = self.xmlRoot.find("reference")
        if refElem is None or "path" not in refElem.attrib:
            return None
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
    
    def getOutputDir(self):
        return self.xmlRoot.attrib["outputDir"]
    
    def setOutputDir(self, path):
        assert os.path.isdir(path)
        self.xmlRoot.attrib["outputDir"] = path
        
    def getHALPath(self):
        halElem = self.xmlRoot.find("hal")
        if halElem is None or "halPath" not in halElem.attrib:
            return None
        return halElem.attrib["halPath"]
    
    def setHALPath(self, path):
        halElem = self.xmlRoot.find("hal")
        if halElem is None:
            halElem = ET.Element("hal")
            self.xmlRoot.append(halElem)
        halElem.attrib["halPath"] = path
        
    def getHALFastaPath(self):
        halElem = self.xmlRoot.find("hal")
        if halElem is None or "fastaPath" not in halElem.attrib:
            return None
        return halElem.attrib["fastaPath"]
    
    def setHALFastaPath(self, path):
        halElem = self.xmlRoot.find("hal")
        if halElem is None:
            halElem = ET.Element("hal")
            self.xmlRoot.append(halElem)
        halElem.attrib["fastaPath"] = path
        
    def setConstraintsFilePath(self, path):
        self.xmlRoot.attrib["constraints"] = path
    
    def getConstraintsFilePath(self):
        if "constraints" not in self.xmlRoot.attrib:
            return None
        return self.xmlRoot.attrib["constraints"]
        
    def getOutputSequenceDir(self):
        if "outputSequenceDir" not in self.xmlRoot.attrib:
            return os.path.join(self.getOutputDir(), "processedSequences")
        return self.xmlRoot.attrib["outputSequenceDir"]
    
    def setOutputSequenceDir(self, path):
        assert os.path.isdir(path)
        self.xmlRoot.attrib["outputSequenceDir"] = path
            
    def getOutgroupEvents(self):
        if self.xmlRoot.attrib.has_key("outgroup_events"):
            return self.xmlRoot.attrib["outgroup_events"].split()
        return []
    
    def setOutgroupEvents(self, outgroupEvents):
        self.xmlRoot.attrib["outgroup_events"] = " ".join(outgroupEvents)
        
    def getConfigPath(self):
        config = self.xmlRoot.attrib["config"]
        if config == 'default':
            config = os.path.join(cactusRootPath(), "cactus_config.xml")
        if config == 'defaultProgressive':
            config = os.path.join(cactusRootPath(), "cactus_progressive_config.xml")
        return config
    
    def setConfigPath(self, path):
        self.xmlRoot.attrib["config"] = path
        
    def getDiskDatabaseString(self):
        conf = self.xmlRoot.find("cactus_disk").find("st_kv_database_conf")
        return ET.tostring(conf).replace("\n", "").replace("\t","")

    # kind of hack to incorporate secondary db interface without
    # revamping anything else.  
    def getSecondaryDBElem(self):
        secondaryConf = self.diskElem.find("secondary_conf")
        if secondaryConf is None:
            return None
        return DbElemWrapper(secondaryConf)

    # kind of hack to incorporate secondary db interface without
    # revamping anything else.  
    def setSecondaryDBElem(self, dbElemWrapper):
        secondaryConf = self.diskElem.find("secondary_conf")
        if secondaryConf is not None:
            self.diskElem.remove(secondaryConf)
        confAttrib = copy.deepcopy(dbElemWrapper.confElem.attrib)
        secondaryConf = ET.SubElement(self.diskElem, "secondary_conf",
                                      confAttrib)
        dbAttrib = copy.deepcopy(dbElemWrapper.getDbElem().attrib)
        secondaryDb = ET.SubElement(secondaryConf, dbElemWrapper.getDbType(),
                                    dbAttrib)
    
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
    
    # Adds an additional single outgroup sequence to the experiment. 
    def addOutgroupSequence(self, ogName, ogDist, ogPath):
        assert ogName is not None
        tree = MultiCactusTree(self.getTree())
        tree.addOutgroup(ogName, ogDist)
        self.xmlRoot.attrib["species_tree"] = NXNewick().writeString(tree)   
        seqs = "%s %s"  % (self.xmlRoot.attrib["sequences"], ogPath)
        self.xmlRoot.attrib["sequences"] = seqs
        self.seqMap[ogName] = ogPath
        self.setOutgroupEvents(self.getOutgroupEvents() + [ogName])

    # return internal structure that maps event names to paths
    def getSequenceMap(self):
        return self.seqMap
