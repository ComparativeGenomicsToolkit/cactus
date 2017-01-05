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
        if "database_name" not in self.dbElem.attrib:
            return None
        return self.dbElem.attrib["database_name"]

    def setDbName(self, name):
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

    @staticmethod
    def createExperimentWrapper(sequences, newickTreeString, outputDir,
                 outgroupGenomes=None,
                 databaseConf=None, configFile=None,
                 halFile=None, fastaFile=None,
                 constraints=None, progressive=False,
                 outputSequenceDir=None):
        #Establish the basics
        rootElem =  ET.Element("cactus_workflow_experiment")
        rootElem.attrib['species_tree'] = newickTreeString
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
        #Outgroup genomes
        if outgroupGenomes != None:
            self.setOutgroupGenomes(outgroupGenomes)
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

    def getTree(self, onlyThisSubtree=False):
        treeString = self.xmlRoot.attrib["species_tree"]
        ret = NXNewick().parseString(treeString, addImpliedRoots = False)
        if onlyThisSubtree:
            # Get a subtree containing only the reference node and its
            # children, rather than a species tree including the
            # outgroups as well
            multiCactus = MultiCactusTree(ret)
            multiCactus.nameUnlabeledInternalNodes()
            multiCactus.computeSubtreeRoots()
            ret = multiCactus.extractSubTree(self.getReferenceNameFromConfig())
        return ret

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

    def getOutgroupGenomes(self):
        if self.xmlRoot.attrib.has_key("outgroup_genomes"):
            return self.xmlRoot.attrib["outgroup_genomes"].split()
        return []

    def setOutgroupGenomes(self, outgroupGenomes):
        self.xmlRoot.attrib["outgroup_genomes"] = " ".join(outgroupGenomes)

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

    def getSequencePath(self, genome):
        """
        Gets the path to the sequence for this genome.
        Returns None if the genome is not found or there is no sequence.
        """
        genomeNodes = self.xmlRoot.findall("genome[@name='%s']" % genome)

        if len(genomeNodes) == 0:
            return None

        assert len(genomeNodes) == 1
        genomeNode = genomeNodes[0]
        return genomeNode.attrib['sequence']

    def setSequencePath(self, genome, seqPath):
        """
        Sets the path to the sequence for this genome.
        """
        genomeNodes = self.xmlRoot.findall("genome[@name='%s']" % genome)

        if len(genomeNodes) == 0:
            # Need to create a new genome element for this sequence
            genomeNode = ET.SubElement(self.xmlRoot, 'genome')
            genomeNode.attrib['name'] = genome
        else:
            assert len(genomeNodes) == 1
            genomeNode = genomeNodes[0]
        genomeNode.attrib['sequence'] = seqPath

    def getGenomesWithSequence(self):
        """
        Return a list of names of genomes in the problem which have sequence.
        """
        genomeNodes = self.xmlRoot.findall("genome")
        return [node.attrib['name'] for node in genomeNodes if 'sequence' in node.attrib]

    def getSequences(self):
        """
        Convenience method for returning the paths to all sequences in the problem
        (in an arbitrary order).
        """
        return [self.getSequencePath(genome) for genome in self.getGenomesWithSequence()]

    def setTree(self, tree):
        """
        Load a new tree.
        """
        # Write the new string to the XML
        treeString = NXNewick().writeString(tree)
        self.xmlRoot.attrib["species_tree"] = treeString

        # Ensure the changes are reflected in the genome elements
        # (adding and deleting elements as necessary).
        genomesInTree = set(tree.getName(id) for id in tree.postOrderTraversal())
        genomeNodes = self.xmlRoot.findall('genome')
        genomeNamesInXML = set(node.attrib['name'] for node in genomeNodes)
        for node in genomeNodes:
            if node.attrib['name'] not in genomesInTree:
                self.xmlRoot.remove(node)

        for genome in genomesInTree:
            if genome not in genomeNamesInXML:
                node = ET.SubElement(self.xmlRoot, 'genome')
                node.attrib['name'] = genome
