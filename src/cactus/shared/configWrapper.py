#!/usr/bin/env python3

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Interface to the cactus config xml file used
to read the progressive-related fields. it's all
considered optional, with default stored as static
members of the configwrapper class

"""
import xml.etree.ElementTree as ET
from xml.dom import minidom
import sys
from cactus.shared.common import findRequiredNode
from cactus.shared.common import getOptionalAttrib

class ConfigWrapper:
    defaultOutgroupStrategy = 'none'
    defaultDoSelf = 'false'
    defaultCoverageFraction = 0
    defaultSingleCopyStrategy = 'none'
    defaultInternalNodePrefix = 'Anc'
    defaultOutgroupThreshold = None
    defaultOutgroupAncestorQualityFraction = 0.75
    defaultMaxParallelSubtrees = 3
    defaultMaxNumOutgroups = 1

    def __init__(self, xmlRoot):
        self.xmlRoot = xmlRoot

    def writeXML(self, path):
        xmlFile = open(path, "w")
        xmlString = ET.tostring(self.xmlRoot, encoding='unicode')
        xmlString = xmlString.replace("\n", "")
        xmlString = xmlString.replace("\t", "")
        xmlString = minidom.parseString(xmlString).toprettyxml()
        xmlFile.write(xmlString)
        xmlFile.close()

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
            strategy == "greedyLeaves" or strategy == "greedyPreference" or \
            strategy == "dynamic"
        return strategy

    def getOutgroupThreshold(self):
        ogElem = self.getOutgroupElem()
        threshold = self.defaultOutgroupThreshold
        if (ogElem is not None and\
            "threshold" in ogElem.attrib and\
            ogElem.attrib["threshold"].lower() != 'none'):
                threshold = int(ogElem.attrib["threshold"])
        return threshold

    def getOutgroupAncestorQualityFraction(self):
        ogElem = self.getOutgroupElem()
        fraction = self.defaultOutgroupAncestorQualityFraction
        if (ogElem is not None and\
            "ancestor_quality_fraction" in ogElem.attrib and\
            ogElem.attrib["ancestor_quality_fraction"].lower() != 'none'):
            fraction = float(ogElem.attrib["ancestor_quality_fraction"])
        return fraction

    def getMaxNumOutgroups(self):
        ogElem = self.getOutgroupElem()
        maxNumOutgroups = self.defaultMaxNumOutgroups
        if (ogElem is not None and\
            "strategy" in ogElem.attrib and\
            "max_num_outgroups" in ogElem.attrib):
            maxNumOutgroups = int(ogElem.attrib["max_num_outgroups"])
        return maxNumOutgroups

    def getDoTrimStrategy(self):
        trimBlastNode = findRequiredNode(self.xmlRoot, "trimBlast")
        if "doTrimStrategy" in trimBlastNode.attrib:
            return trimBlastNode.attrib["doTrimStrategy"] == "1"
        return False

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

    def getKtserverMemory(self, default=sys.maxsize):
        ktServerElem = self.xmlRoot.find("ktserver")
        if ktServerElem is not None and "memory" in ktServerElem.attrib:
            return int(ktServerElem.attrib["memory"])
        return default

    def getKtserverCpu(self, default=sys.maxsize):
        ktServerElem = self.xmlRoot.find("ktserver")
        if ktServerElem is not None and "cpu" in ktServerElem.attrib:
            return int(ktServerElem.attrib["cpu"])
        return default

    def getDefaultMemory(self):
        constantsElem = self.xmlRoot.find("constants")
        return int(constantsElem.attrib["defaultMemory"])

    def getExportHalDisk(self):
        exportHalElem = self.xmlRoot.find("exportHal")
        return int(exportHalElem.attrib["disk"])

    def substituteAllPredefinedConstantsWithLiterals(self):
        constants = findRequiredNode(self.xmlRoot, "constants")
        defines = constants.find("defines")
        def replaceAllConstants(node, defines):
            for attrib in node.attrib:
                if node.attrib[attrib] in defines.attrib:
                    node.attrib[attrib] = defines.attrib[node.attrib[attrib]]
            for child in node:
                replaceAllConstants(child, defines)
        if defines != None:
            replaceAllConstants(self.xmlRoot, defines)
            constants.remove(defines)

    def substituteAllDivergenceContolledParametersWithLiterals(self, maxDivergence):
        constants = findRequiredNode(self.xmlRoot, "constants")
        divergences = constants.find("divergences")
        messages = []
        if divergences != None:
            useDefaultDivergences = getOptionalAttrib(divergences, attribName="useDefault", typeFn=bool, default=False)
            def replaceAllDivergenceParameters(node):
                for child in node:
                    if child.tag == "divergence":
                        attribName = child.attrib["argName"]
                        arg = child.attrib["default"]
                        divergence = sys.maxsize
                        if not useDefaultDivergences:
                            for i in list(child.attrib.keys()):
                                if i in list(divergences.attrib.keys()):
                                    j = float(divergences.attrib[i])
                                    if j < divergence and j >= maxDivergence:
                                        arg = child.attrib[i]
                                        divergence = j
                        messages.append("Made argument %s=%s in tag %s with divergence threshold of %s for longest path of %s (useDefaultDivergences=%s)" % (attribName, arg, node.tag, divergence, maxDivergence, useDefaultDivergences))
                        node.attrib[attribName] = arg
                    else:
                        replaceAllDivergenceParameters(child)
            replaceAllDivergenceParameters(self.xmlRoot)
        return messages

    def turnAllModesOn(self):
        """Switches on check, normalisation etc. to use when debugging/testing
        """
        findRequiredNode(self.xmlRoot, "check").attrib["runCheck"] = "1"
        findRequiredNode(self.xmlRoot, "normal").attrib["iterations"] = "2"

    def turnOffHeaderChecks(self):
        """Turns off the preprocessor stage that checks whether headers can be
        used in an assembly hub."""
        for node in self.xmlRoot.findall("preprocessor"):
            if 'checkAssemblyHub' in node.attrib:
                node.attrib['checkAssemblyHub'] = '0'
