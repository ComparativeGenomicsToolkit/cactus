#!/usr/bin/env python3

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Basic interface to the multi cactus project xml file.

"""
import xml.etree.ElementTree as ET
from xml.dom import minidom

from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.experimentWrapper import ExperimentWrapper
from sonLib.nxnewick import NXNewick

from toil.lib.bioio import logger

class MultiCactusProject:
    def __init__(self):
        self.mcTree = None
        self.expMap = dict()
        self.expIDMap = None
        self.inputSequenceMap = {}
        self.inputSequenceIDMap = {}
        self.outputSequenceIDMap = {}
        self.configID = None

    def readXML(self, path):
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
                self.expIDMap[nameElem] = cactusPathElem.attrib["experiment_id"]
        self.inputSequenceMap = dict(list(zip(xmlRoot.attrib["inputSequenceNames"].split(),
                                         xmlRoot.attrib["inputSequences"].split())))
        if "inputSequenceIDs" in xmlRoot.attrib:
            self.inputSequenceIDMap = dict(list(zip(xmlRoot.attrib["inputSequenceIDNames"].split(),
                                               xmlRoot.attrib["inputSequenceIDs"].split())))
        if "outputSequenceIDs" in xmlRoot.attrib:
            self.outputSequenceIDMap = dict(list(zip(xmlRoot.attrib["outputSequenceNames"].split(),
                                                xmlRoot.attrib["outputSequenceIDs"].split())))

        logger.info("xmlRoot = %s" % ET.tostring(xmlRoot, encoding='unicode'))
        if "configID" in xmlRoot.attrib:
            self.configID = xmlRoot.attrib["configID"]

        self.mcTree.assignSubtreeRootNames(self.expMap)

    def writeXML(self, path):
        xmlRoot = ET.Element("multi_cactus")
        treeElem = ET.Element("tree")
        treeElem.text = NXNewick().writeString(self.mcTree)
        xmlRoot.append(treeElem)
        for name, expPath in list(self.expMap.items()):
            cactusPathElem = ET.Element("cactus")
            cactusPathElem.attrib["name"] = name
            cactusPathElem.attrib["experiment_path"] = expPath
            if self.expIDMap:
                cactusPathElem.attrib["experiment_id"] = self.expIDMap[name]
            xmlRoot.append(cactusPathElem)
        #We keep track of all the input sequences at the top level
        xmlRoot.attrib["inputSequences"] = " ".join(list(self.inputSequenceMap.values()))
        xmlRoot.attrib["inputSequenceNames"] = " ".join(list(self.inputSequenceMap.keys()))
        if self.inputSequenceIDMap:
            xmlRoot.attrib["inputSequenceIDs"] = " ".join(list(self.inputSequenceIDMap.values()))
            xmlRoot.attrib["inputSequenceIDNames"] = " ".join(list(self.inputSequenceIDMap.keys()))
        if self.outputSequenceIDMap:
            xmlRoot.attrib["outputSequenceIDs"] = " ".join(list(self.outputSequenceIDMap.values()))
            xmlRoot.attrib["outputSequenceNames"] = " ".join(list(self.outputSequenceIDMap.keys()))
        if self.configID:
            xmlRoot.attrib["configID"] = self.configID

        xmlFile = open(path, "w")
        xmlString = ET.tostring(xmlRoot, encoding='unicode')
        xmlString = minidom.parseString(xmlString).toprettyxml()
        xmlFile.write(xmlString)
        xmlFile.close()

    def syncToFileStore(self, toil):
        self.expIDMap = dict()
        for name, expPath in list(self.expMap.items()):
            expWrapper = ExperimentWrapper(ET.parse(expPath).getroot())
            expWrapper.setConfigID(toil.importFile("file://" + expWrapper.getConfigPath()))
            expWrapper.writeXML(expPath)
            self.expIDMap[name] = toil.importFile("file://" + expPath)

    def getConfigPath(self):
        return ExperimentWrapper(ET.parse(list(self.expMap.values())[0]).getroot()).getConfigPath()

    def setConfigID(self, configID):
        self.configID = configID

    def getConfigID(self):
        return self.configID

if __name__ == '__main__':
    main()
