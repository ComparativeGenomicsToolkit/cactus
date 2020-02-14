#!/usr/bin/env python3

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

"""Interface to the cactus experiment xml file used
to read and modify an existing experiment"""
import os
import xml.etree.ElementTree as ET
from xml.dom import minidom

from cactus.progressive.multiCactusTree import MultiCactusTree
from sonLib.nxnewick import NXNewick
from cactus.shared.common import cactusRootPath

class DbElemWrapper(object):
    def __init__(self, confElem):
        typeString = confElem.attrib["type"]
        dbElem = confElem.find(typeString)
        self.dbElem = dbElem
        self.confElem = confElem
        self.dbElem.attrib["database_dir"] = "fakepath"

    def check(self):
        """Function checks the database conf is as expected and creates useful exceptions
        if not"""
        dataString = ET.tostring(self.confElem, encoding='unicode')
        if self.confElem.tag != "st_kv_database_conf":
            raise RuntimeError("The database conf string is improperly formatted: %s" % dataString)
        if "type" not in self.confElem.attrib:
            raise RuntimeError("The database conf string does not have a type attrib: %s" % dataString)
        typeString = self.confElem.attrib["type"]
        if typeString == "tokyo_cabinet":
            tokyoCabinet = self.confElem.find("tokyo_cabinet")
            if tokyoCabinet == None:
                raise RuntimeError("Database conf is of type tokyo cabinet but there is no nested tokyo cabinet tag: %s" % dataString)
            if "database_dir" not in tokyoCabinet.attrib:
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
        return ET.tostring(self.confElem, encoding='unicode')

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

class ExperimentWrapper(DbElemWrapper):
    def __init__(self, xmlRoot):
        self.diskElem = xmlRoot.find("cactus_disk")
        confElem = self.diskElem.find("st_kv_database_conf")
        super(ExperimentWrapper, self).__init__(confElem)
        self.xmlRoot = xmlRoot

    @staticmethod
    def createExperimentWrapper(newickTreeString,
                                genomes,
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
            ET.SubElement(databaseConf, "tokyo_cabinet")
        self = ExperimentWrapper(rootElem)
        #Setup the config
        self.setConfigPath("default")
        if progressive == True:
            self.setConfigPath("defaultProgressive")
        if configFile != None:
            self.setConfigPath(configFile)
        for genome in genomes:
            genomeNode = ET.SubElement(rootElem, "genome")
            genomeNode.attrib['name'] = genome
        #Constraints
        if constraints != None:
            self.setConstraintsFilePath(constraints)
        #Outgroup genomes
        if outgroupGenomes != None:
            self.setOutgroupGenomes(outgroupGenomes)
        return self

    def writeXML(self, path):
        #Replacement for writeExperimentFile
        xmlFile = open(path, "w")
        xmlString = ET.tostring(self.xmlRoot, encoding='unicode')
        xmlString = xmlString.replace("\n", "")
        xmlString = xmlString.replace("\t", "")
        xmlString = minidom.parseString(xmlString).toprettyxml()
        xmlFile.write(xmlString)
        xmlFile.close()

    def setConfigID(self, configID):
        self.xmlRoot.attrib["configID"] = str(configID)

    def getConfigID(self):
        return self.xmlRoot.attrib["configID"]

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
            ret = multiCactus.extractSubTree(self.getRootGenome())
        return ret

    def isRootReconstructed(self):
        """
        Return True if this is a reconstruction problem and False otherwise.
        """
        rootElem = self.xmlRoot.find('root')
        if rootElem is None:
            return False
        return 'reconstruction' in rootElem.attrib and rootElem.attrib['reconstruction'] == '1'

    def setRootReconstructed(self, reconstructed):
        """
        Set whether this is a reconstruction problem or not.
        """
        rootElem = self.xmlRoot.find('root')
        if reconstructed:
            rootElem.attrib['reconstruction'] = '1'
        else:
            refElem = self.xmlRoot.find("reference")
            if refElem is not None:
                del refElem
            if 'reconstruction' in rootElem.attrib:
                del rootElem.attrib['reconstruction']

    def getRootGenome(self):
        """
        Get the root of the subtree to be aligned.
        """
        rootElem = self.xmlRoot.find('root')
        if rootElem is None:
            return None
        # The text can contain '\n's and whitespace so we remove them
        return rootElem.text.strip()

    def setRootGenome(self, root):
        """
        Set the root of the subtree to be aligned.

        (the genome to be reconstructed if this is a reconstruction problem)
        """
        rootElem = self.xmlRoot.find('root')
        if rootElem is None:
            rootElem = ET.SubElement(self.xmlRoot, 'root')
        rootElem.text = root

    def setReferenceID(self, refID):
        '''Set the file store ID of the reconstructed ancestral
        genome for this experiment. This should be downloaded
        onto the master node after the experiment has finished running.'''
        refElem = self.xmlRoot.find("reference")
        if refElem is None:
            refElem = ET.SubElement(self.xmlRoot, "reference")
        refElem.attrib["id"] = str(refID)

    def getReferenceID(self):
        refElem = self.xmlRoot.find("reference")
        if refElem is not None and "id" in refElem.attrib:
            return refElem.attrib["id"]
        else:
            return None

    def setHalID(self, halID):
        '''Set the file store ID of the HAL file
        resulting from this experiment.'''
        halElem = self.xmlRoot.find("hal")
        if halElem is None:
            halElem = ET.SubElement(self.xmlRoot, "hal")
        halElem.attrib["halID"] = str(halID)

    def getHalID(self):
        halElem = self.xmlRoot.find("hal")
        return halElem.attrib["halID"]

    def setHalFastaID(self, halFastaID):
        halElem = self.xmlRoot.find("hal")
        if halElem is None:
            halElem = ET.SubElement(self.xmlRoot, "hal")
        halElem.attrib["fastaID"] = str(halFastaID)

    def getHalFastaID(self):
        halElem = self.xmlRoot.find("hal")
        return halElem.attrib["fastaID"]

    def getOutgroupGenomes(self):
        genomeNodes = self.xmlRoot.findall("genome")
        return [genome.attrib['name'] for genome in genomeNodes if 'outgroup' in genome.attrib and genome.attrib['outgroup'] == "1"]

    def setOutgroupGenomes(self, outgroupGenomes):
        genomeNodes = self.xmlRoot.findall("genome")
        for node in genomeNodes:
            if 'outgroup' in node.attrib:
                del node.attrib['outgroup']
        outgroupNodes = [node for node in genomeNodes if node.attrib['name'] in outgroupGenomes]
        for node in outgroupNodes:
            node.attrib['outgroup'] = "1"

    def setConstraintsID(self, fileID):
        self.xmlRoot.attrib["constraintsID"] = str(fileID)

    def getConstraintsID(self, fileID):
        return self.xmlRoot.attrib["constraintsID"]

    def getConfigPath(self):
        config = self.xmlRoot.attrib["config"]
        if config == 'default':
            config = os.path.join(cactusRootPath(), "cactus_config.xml")
        if config == 'defaultProgressive':
            config = os.path.join(cactusRootPath(), "cactus_progressive_config.xml")
        return config

    def setConfigPath(self, path):
        self.xmlRoot.attrib["config"] = path

    def getSequenceID(self, genome):
        """
        Gets the fileStoreID for the sequence for this genome.
        Returns None if the genome is not found or there is no sequence.
        """
        genomeNodes = self.xmlRoot.findall("genome[@name='%s']" % genome)

        if len(genomeNodes) == 0:
            return None

        assert len(genomeNodes) == 1
        genomeNode = genomeNodes[0]
        if 'sequence' in genomeNode.attrib:
            id = genomeNode.attrib['sequence']
        else:
            id = None
        return id

    def setSequenceID(self, genome, seqID):
        """
        Sets the fileStoreID for the sequence for this genome.
        """
        seqID = str(seqID)
        genomeNodes = self.xmlRoot.findall("genome[@name='%s']" % genome)

        if len(genomeNodes) == 0:
            # Need to create a new genome element for this sequence
            genomeNode = ET.SubElement(self.xmlRoot, 'genome')
            genomeNode.attrib['name'] = genome
        else:
            assert len(genomeNodes) == 1
            genomeNode = genomeNodes[0]
        genomeNode.attrib['sequence'] = seqID

    def getGenomesWithSequence(self):
        """
        Return a list of names of genomes in the problem which have sequence.
        """
        genomeNodes = self.xmlRoot.findall("genome")
        return [node.attrib['name'] for node in genomeNodes if 'sequence' in node.attrib]

    def getSequenceIDs(self):
        """
        Convenience method for returning the paths to all sequences in the problem
        (in an arbitrary order).
        """
        return [self.getSequenceID(genome) for genome in self.getGenomesWithSequence()]

    def setTree(self, tree):
        """
        Load a new tree.
        """
        # Write the new string to the XML
        treeString = NXNewick().writeString(tree)
        self.xmlRoot.attrib["species_tree"] = treeString

        # Ensure the changes are reflected in the genome elements
        # (adding and deleting elements as necessary).
        genomesInTree = set(tree.getName(id) for id in tree.postOrderTraversal() if tree.hasName(id))
        genomeNodes = self.xmlRoot.findall('genome')
        genomeNamesInXML = set(node.attrib['name'] for node in genomeNodes)
        for node in genomeNodes:
            if node.attrib['name'] not in genomesInTree:
                self.xmlRoot.remove(node)

        for genome in genomesInTree:
            if genome not in genomeNamesInXML:
                node = ET.SubElement(self.xmlRoot, 'genome')
                node.attrib['name'] = genome
