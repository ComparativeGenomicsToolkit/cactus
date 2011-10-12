#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Interface to the cactus config xml file used
to read the progressive-related fields. it's all 
considered optional, with default stored as static
members of the configwrapper class

"""
import unittest

import os
import xml.etree.ElementTree as ET
from xml.dom import minidom
import sys
import random
import math
import copy


class ConfigWrapper:
    defaultOutgroupStrategy = 'none'
    defaultSubtreeSize = 2
    defaultDoSelf = 'false'
    defaultCoverageFraction = 0
    defaultSingleCopyStrategy = 'none'
    defaultInternalNodePrefix = 'Anc'
    
    def __init__(self, xmlRoot):
        self.xmlRoot = xmlRoot

    def writeXML(self, path):
        xmlFile = open(path, "w")
        xmlString = ET.tostring(self.xmlRoot)
        xmlString = xmlString.replace("\n", "")
        xmlString = xmlString.replace("\t", "")
        xmlString = minidom.parseString(xmlString).toprettyxml()
        xmlFile.write(xmlString)
        xmlFile.close()
    
    def setReferenceName(self, name):
        refElem = self.xmlRoot.find("reference")
        refElem.attrib["reference"] = name
        
    def getMCElem(self):
        return self.xmlRoot.find("multi_cactus")
    
    def getOutgroupElem(self):
        mcElem = self.getMCElem()
        if mcElem is not None:
            return mcElem.find("outgroup")
    
    def getCoverageElem(self):
        mcElem = self.getMCElem()
        if mcElem is not None:
            return mcElem.find("coverage")
    
    def getDecompositionElem(self):
        mcElem = self.getMCElem()
        if mcElem is not None:
            return mcElem.find("decomposition")
    
    def getOutgroupStrategy(self):
        ogElem = self.getOutgroupElem()
        strategy = self.defaultOutgroupStrategy
        if ogElem is not None and "strategy" in ogElem.attrib:
            strategy = ogElem.attrib["strategy"]
        assert strategy == "none" or strategy == "greedy" or \
            strategy == "greedyLeaves"
        return strategy
    
    def getSubtreeSize(self):
        decompElem = self.getDecompositionElem()
        subtreeSize = self.defaultSubtreeSize
        if decompElem is not None and "subtree_size" in decompElem.attrib:
            subtreeSize = int(decompElem.attrib["subtree_size"])
        assert subtreeSize > 1
        return subtreeSize
            
    def getDoSelfAlignment(self):
        decompElem = self.getDecompositionElem()
        doSelf = self.defaultDoSelf
        if decompElem is not None and "self_alignment" in decompElem.attrib:
            doSelf = decompElem.attrib["self_alignment"].lower()
        assert doSelf == "true" or doSelf == "false"
        return doSelf == "true"
    
    def getDefaultInternalNodePrefix(self):
        decompElem = self.getDecompositionElem()
        prefix = self.defaultInternalNodePrefix
        if decompElem is not None and\
         "default_internal_node_prefix" in decompElem.attrib:
            prefix = decompElem.attrib["default_internal_node_prefix"]
        assert len(prefix) > 0
        return prefix
    
    def getRequiredFraction(self):
        coverageElem = self.getCoverageElem()
        frac = self.defaultCoverageFraction
        if coverageElem is not None \
        and "required_fraction" in coverageElem.attrib:
            frac = float(coverageElem.attrib["required_fraction"])
        assert frac >= 0 and frac <= 1
        return frac
    
    def getSingleCopyStrategy(self):
        coverageElem = self.getCoverageElem()
        singleCopy = self.defaultSingleCopyStrategy
        if coverageElem is not None \
        and "single_copy_strategy" in coverageElem.attrib:
            singleCopy = coverageElem.attrib["single_copy_strategy"]
        assert singleCopy == 'none' or singleCopy == 'outgroup' \
        or singleCopy == 'all'
        return singleCopy
        
        