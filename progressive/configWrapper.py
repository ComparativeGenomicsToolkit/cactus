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
    defaultOutgroupThreshold = None
    defaultOutgroupAncestorQualityFraction = 0.75
    defaultMaxParallelSubtrees = 3
    
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
    
    def getOutgroupThreshold(self):
        ogElem = self.getOutgroupElem()
        threshold = self.defaultOutgroupThreshold
        if (ogElem is not None and\
            "strategy" in ogElem.attrib and\
            ogElem.attrib["strategy"] == "greedy" and\
            "threshold" in ogElem.attrib and\
            ogElem.attrib["threshold"].lower() != 'none'):
                threshold = int(ogElem.attrib["threshold"])
        return threshold

    def getOutgroupAncestorQualityFraction(self):
        ogElem = self.getOutgroupElem()
        fraction = self.defaultOutgroupAncestorQualityFraction
        if (ogElem is not None and\
            "strategy" in ogElem.attrib and\
            ogElem.attrib["strategy"] == "greedy" and\
            "ancestor_quality_fraction" in ogElem.attrib and\
            ogElem.attrib["ancestor_quality_fraction"].lower() != 'none'):
            fraction = float(ogElem.attrib["ancestor_quality_fraction"])
        return fraction
    
    def getSubtreeSize(self):
        decompElem = self.getDecompositionElem()
        subtreeSize = self.defaultSubtreeSize
        if decompElem is not None and "subtree_size" in decompElem.attrib:
            subtreeSize = int(decompElem.attrib["subtree_size"])
        assert subtreeSize > 1
        return subtreeSize

    def setSubtreeSize(self, subtreeSize):
        decompElem = self.getDecompositionElem()
        assert decompElem is not None
        decompElem.attrib["subtree_size"] = str(subtreeSize)
            
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
    
    def getBuildHal(self):
        halElem = self.xmlRoot.find("hal")
        if halElem is not None and "buildHal" in halElem.attrib:
            build = halElem.attrib["buildHal"]
            if build == "1" or build.lower() == "true":
                return True
        return False

    def setBuildHal(self, buildHal):
         halElem = self.xmlRoot.find("hal")
         assert halElem is not None
         halElem.attrib["buildHal"] = str(int(buildHal))

    def getBuildMaf(self):
        halElem = self.xmlRoot.find("hal")
        if halElem is not None and "buildMaf" in halElem.attrib:
            build = halElem.attrib["buildMaf"]
            if build == "1" or build.lower() == "true":
                return True
        return False

    def setBuildMaf(self, buildMaf):
         halElem = self.xmlRoot.find("hal")
         assert halElem is not None
         halElem.attrib["buildMaf"] = str(int(buildMaf))

    def getBuildFasta(self):
        halElem = self.xmlRoot.find("hal")
        if halElem is not None and "buildFasta" in halElem.attrib:
            build = halElem.attrib["buildFasta"]
            if build == "1" or build.lower() == "true":
                return True
        return False

    def setBuildFasta(self, buildFasta):
        halElem = self.xmlRoot.find("hal")
        assert halElem is not None
        halElem.attrib["buildFasta"] = str(int(buildFasta))
    
    def getJoinMaf(self):
        halElem = self.xmlRoot.find("hal")
        if halElem is not None and "joinMaf" in halElem.attrib:            
            maf = halElem.attrib["joinMaf"]
            if maf == "0" or maf.lower() == "false":
                return False
        return self.getBuildMaf()

    def setJoinMaf(self, joinMaf):
         halElem = self.xmlRoot.find("hal")
         assert halElem is not None
         halElem.attrib["joinMaf"] = str(int(joinMaf))

    def getMaxParallelSubtrees(self):
        decompElem = self.getDecompositionElem()
        maxParallelSubtrees = self.defaultMaxParallelSubtrees
        if decompElem is not None and\
               "max_parallel_subtrees" in decompElem.attrib:
            maxParallelSubtrees = int(
                decompElem.attrib["max_parallel_subtrees"])
        assert maxParallelSubtrees > 0
        return maxParallelSubtrees

    def setMaxParallelSubtrees(self, maxParallel):
        decompElem = self.getDecompositionElem()
        assert decompElem is not None
        decompElem.attrib["max_parallel_subtrees"] = str(maxParallel)
            
    # the minBlockDegree, when specified in the final, "base" 
    # iteration, does not play nicely with the required fraction
    # option
    def verifyMinBlockDegree(self, experiment):  
        barElem = self.xmlRoot.find("bar")
        minBlockDegree = int(barElem.attrib["minimumBlockDegree"])
        reqSpecies = experiment.xmlRoot.find("required_species")
        if reqSpecies is not None and minBlockDegree != 2:
            sys.stderr.write("WARNING: reverting minBlockDegree from %i to 2 because of required_species\n" % 
                             minBlockDegree)
            itElem.attrib["minBlockDegree"] = "2"
            
        
