#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Script strings together all the components to make the basic pipeline for reconstruction.

The script uses the the jobTree.scriptTree target framework to structure all the related wrappers.
"""

import os
import sys
import xml.etree.ElementTree as ET
import math
import time
import bz2
import random
import copy
import shutil
import networkx as NX
from optparse import OptionParser

from sonLib.bioio import getTempFile
from sonLib.bioio import newickTreeParser

from sonLib.bioio import logger
from sonLib.bioio import setLoggingFromOptions
from sonLib.bioio import system, popenCatch
from sonLib.bioio import makeSubDir
from sonLib.nxnewick import NXNewick

from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.common import cactusRootPath
  
from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import getLogLevelString

from cactus.shared.common import getOptionalAttrib
from cactus.shared.common import runCactusSetup
from cactus.shared.common import runCactusCaf
from cactus.shared.common import runCactusGetFlowers
from cactus.shared.common import runCactusExtendFlowers
from cactus.shared.common import runCactusSplitFlowersBySecondaryGrouping
from cactus.shared.common import encodeFlowerNames
from cactus.shared.common import decodeFirstFlowerName
from cactus.shared.common import runCactusConvertAlignmentToCactus
from cactus.shared.common import runCactusPhylogeny
from cactus.shared.common import runCactusAdjacencies
from cactus.shared.common import runCactusBar
from cactus.shared.common import runCactusMakeNormal 
from cactus.shared.common import runCactusReference
from cactus.shared.common import runCactusAddReferenceCoordinates
from cactus.shared.common import runCactusCheck
from cactus.shared.common import runCactusHalGenerator
from cactus.shared.common import runCactusFlowerStats
from cactus.shared.common import runCactusSecondaryDatabase
from cactus.shared.common import runCactusFastaGenerator
from cactus.shared.common import findRequiredNode
from cactus.shared.common import runConvertAlignmentsToInternalNames
from cactus.shared.common import runStripUniqueIDs

from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.blast.cactus_blast import BlastIngroupsAndOutgroups
from cactus.blast.cactus_blast import BlastFlower
from cactus.blast.cactus_blast import BlastOptions

from cactus.preprocessor.cactus_preprocessor import CactusPreprocessor

from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.shared.experimentWrapper import DbElemWrapper
from cactus.shared.configWrapper import ConfigWrapper
from cactus.pipeline.ktserverJobTree import addKtserverDependentChild

############################################################
############################################################
############################################################
##Shared functions
############################################################
############################################################
############################################################

def extractNode(node):
    """Make an XML node free of its parent subtree
    """
    return ET.fromstring(ET.tostring(node))

def getTargetNode(phaseNode, targetClass):
    """Gets a target node for a given target.
    """
    className = targetClass.__name__
    assert className != ''
    assert className.isalnum()
    return phaseNode.find(className)

class CactusTarget(Target):
    """Base target for all cactus workflow targets.
    """
    def __init__(self, phaseNode, constantsNode, overlarge=False):
        self.phaseNode = phaseNode
        self.constantsNode = constantsNode
        self.overlarge = overlarge
        self.targetNode = getTargetNode(self.phaseNode, self.__class__)
        if overlarge:
            Target.__init__(self, memory=self.getOptionalTargetAttrib("overlargeMemory", typeFn=int, 
                                                                      default=getOptionalAttrib(self.constantsNode, "defaultOverlargeMemory", int, default=sys.maxint)),
                                  cpu=self.getOptionalTargetAttrib("overlargeCpu", typeFn=int, 
                                                                      default=getOptionalAttrib(self.constantsNode, "defaultOverlargeCpu", int, default=sys.maxint)))
        else:
            Target.__init__(self, memory=self.getOptionalTargetAttrib("memory", typeFn=int, 
                                                                      default=getOptionalAttrib(self.constantsNode, "defaultMemory", int, default=sys.maxint)),
                                  cpu=self.getOptionalTargetAttrib("cpu", typeFn=int, 
                                                                      default=getOptionalAttrib(self.constantsNode, "defaultCpu", int, default=sys.maxint)))
    
    def getOptionalPhaseAttrib(self, attribName, typeFn=None, default=None):
        """Gets an optional attribute of the phase node.
        """
        return getOptionalAttrib(node=self.phaseNode, attribName=attribName, typeFn=typeFn, default=default)
    
    def getOptionalTargetAttrib(self, attribName, typeFn=None, default=None):
        """Gets an optional attribute of the target node.
        """
        return getOptionalAttrib(node=self.targetNode, attribName=attribName, typeFn=typeFn, default=default)

class CactusPhasesTarget(CactusTarget):
    """Base target for each workflow phase target.
    """
    def __init__(self, cactusWorkflowArguments, phaseName, topFlowerName=0, index=0):
        phaseNode = findRequiredNode(cactusWorkflowArguments.configNode, phaseName, index)
        constantsNode = findRequiredNode(cactusWorkflowArguments.configNode, "constants")
        CactusTarget.__init__(self, phaseNode=phaseNode, constantsNode=constantsNode, overlarge=False)
        self.index = index
        self.cactusWorkflowArguments = cactusWorkflowArguments
        self.topFlowerName = topFlowerName
    
    def makeRecursiveChildTarget(self, target, launchSecondaryKtForRecursiveTarget=False):
        newChild = target(phaseNode=extractNode(self.phaseNode), 
                          constantsNode=extractNode(self.constantsNode),
                          cactusDiskDatabaseString=self.cactusWorkflowArguments.cactusDiskDatabaseString, 
                          flowerNames=encodeFlowerNames((self.topFlowerName,)), overlarge=True,
                          cactusWorkflowArguments=self.cactusWorkflowArguments)

        if launchSecondaryKtForRecursiveTarget and ExperimentWrapper(self.cactusWorkflowArguments.experimentNode).getDbType() == "kyoto_tycoon":
            cw = ConfigWrapper(self.cactusWorkflowArguments.configNode)
            memory = cw.getKtserverMemory(default=getOptionalAttrib(
                    self.constantsNode, "defaultMemory", int, default=sys.maxint))
            cpu = cw.getKtserverCpu(default=getOptionalAttrib(
                    self.constantsNode, "defaultCpu", int, default=sys.maxint))
            addKtserverDependentChild(self, newChild, maxMemory=memory, maxCpu=cpu, isSecondary = True)
        else:
            self.addChildTarget(newChild)
    
    def makeFollowOnPhaseTarget(self, target, phaseName, index=0):
        self.setFollowOnTarget(target(cactusWorkflowArguments=self.cactusWorkflowArguments, phaseName=phaseName, 
                                      topFlowerName=self.topFlowerName, index=index))
        
    def runPhase(self, recursiveTarget, nextPhaseTarget, nextPhaseName, doRecursion=True, index=0, launchSecondaryKtForRecursiveTarget=False):
        self.logToMaster("Starting %s phase target with index %i at %s seconds (recursing = %i)" % (self.phaseNode.tag, self.getPhaseIndex(), time.time(), doRecursion))
        if doRecursion:
            self.makeRecursiveChildTarget(recursiveTarget, launchSecondaryKtForRecursiveTarget)
        self.makeFollowOnPhaseTarget(target=nextPhaseTarget, phaseName=nextPhaseName, index=index)
        
    def getPhaseIndex(self):
        return self.index
    
    def getPhaseNumber(self):
        return len(self.cactusWorkflowArguments.configNode.findall(self.phaseNode.tag))
    
    def setupSecondaryDatabase(self):
        """Setup the secondary database
        """
        confXML = ET.fromstring(self.cactusWorkflowArguments.secondaryDatabaseString)
        dbElem = DbElemWrapper(confXML)
        if dbElem.getDbType() != "kyoto_tycoon":
            runCactusSecondaryDatabase(self.cactusWorkflowArguments.secondaryDatabaseString, create=True)
    
    def cleanupSecondaryDatabase(self):
        """Cleanup the secondary database
        """
        confXML = ET.fromstring(self.cactusWorkflowArguments.secondaryDatabaseString)
        dbElem = DbElemWrapper(confXML)
        if dbElem.getDbType() != "kyoto_tycoon":
            runCactusSecondaryDatabase(self.cactusWorkflowArguments.secondaryDatabaseString, create=False)

class CactusRecursionTarget(CactusTarget):
    """Base recursive target for traversals up and down the cactus tree.
    """
    maxSequenceSizeOfFlowerGroupingDefault = 1000000
    def __init__(self, phaseNode, constantsNode, cactusDiskDatabaseString, flowerNames, overlarge=False, cactusWorkflowArguments=None):
        CactusTarget.__init__(self, phaseNode=phaseNode, constantsNode=constantsNode, overlarge=overlarge)
        self.cactusDiskDatabaseString = cactusDiskDatabaseString
        self.flowerNames = flowerNames
        self.cactusWorkflowArguments = cactusWorkflowArguments
        
    def makeFollowOnRecursiveTarget(self, target, phaseNode=None):
        """Sets the followon to the given recursive target
        """
        if phaseNode == None:
            phaseNode = self.phaseNode
        self.setFollowOnTarget(target(phaseNode=phaseNode, constantsNode=self.constantsNode,
                                      cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                                      flowerNames=self.flowerNames, overlarge=self.overlarge,
                                      cactusWorkflowArguments=self.cactusWorkflowArguments))
        
    def makeChildTargets(self, flowersAndSizes, target, overlargeTarget=None, 
                         phaseNode=None, runFlowerStats=False):
        """Make a set of child targets for a given set of flowers and chosen child target
        """
        if overlargeTarget == None:
            overlargeTarget = target
        if phaseNode == None:
            phaseNode = self.phaseNode
        for overlarge, flowerNames in flowersAndSizes:
            if overlarge: #Make sure large flowers are on there own, in their own job
                if runFlowerStats:
                    flowerStatsString = runCactusFlowerStats(cactusDiskDatabaseString=self.cactusDiskDatabaseString, flowerName=decodeFirstFlowerName(flowerNames))
                    self.logToMaster("Adding an oversize flower for target class %s and stats %s" \
                                             % (overlargeTarget, flowerStatsString))
                else:
                    self.logToMaster("Adding an oversize flower %s for target class %s" \
                                             % (decodeFirstFlowerName(flowerNames), overlargeTarget))
                self.addChildTarget(overlargeTarget(cactusDiskDatabaseString=self.cactusDiskDatabaseString, phaseNode=phaseNode, 
                                                    constantsNode=self.constantsNode,
                                                    flowerNames=flowerNames, overlarge=True, #This ensures overlarge flowers,
                                                    cactusWorkflowArguments=self.cactusWorkflowArguments))
            else:
                self.addChildTarget(target(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                                           phaseNode=phaseNode, constantsNode=self.constantsNode, flowerNames=flowerNames, overlarge=False,
                                           cactusWorkflowArguments=self.cactusWorkflowArguments))
        
    def makeRecursiveTargets(self, target=None, phaseNode=None, runFlowerStats=False):
        """Make a set of child targets for a given set of parent flowers.
        """
        if target == None:
            target = self.__class__
        targetNode = getTargetNode(self.phaseNode, target)
        flowersAndSizes=runCactusGetFlowers(cactusDiskDatabaseString=self.cactusDiskDatabaseString, flowerNames=self.flowerNames, 
                                            minSequenceSizeOfFlower=getOptionalAttrib(targetNode, "minFlowerSize", int, 0), 
                                            maxSequenceSizeOfFlowerGrouping=getOptionalAttrib(targetNode, "maxFlowerGroupSize", int, 
                                            default=CactusRecursionTarget.maxSequenceSizeOfFlowerGroupingDefault),
                                            maxSequenceSizeOfSecondaryFlowerGrouping=getOptionalAttrib(targetNode, "maxFlowerWrapperGroupSize", int, 
                                            default=CactusRecursionTarget.maxSequenceSizeOfFlowerGroupingDefault))
        self.makeChildTargets(flowersAndSizes=flowersAndSizes, 
                              target=target, phaseNode=phaseNode,
                              runFlowerStats=runFlowerStats)
    
    def makeExtendingTargets(self, target, overlargeTarget=None, phaseNode=None, runFlowerStats=False):
        """Make set of child targets that extend the current cactus tree.
        """
        targetNode = getTargetNode(self.phaseNode, target)
        flowersAndSizes=runCactusExtendFlowers(cactusDiskDatabaseString=self.cactusDiskDatabaseString, flowerNames=self.flowerNames, 
                                              minSequenceSizeOfFlower=getOptionalAttrib(targetNode, "minFlowerSize", int, 0), 
                                              maxSequenceSizeOfFlowerGrouping=getOptionalAttrib(targetNode, "maxFlowerGroupSize", int, 
                                              default=CactusRecursionTarget.maxSequenceSizeOfFlowerGroupingDefault))
        self.makeChildTargets(flowersAndSizes=flowersAndSizes, 
                              target=target, overlargeTarget=overlargeTarget, 
                              phaseNode=phaseNode,
                              runFlowerStats=runFlowerStats)
    
    def makeWrapperTargets(self, target, overlargeTarget=None, phaseNode=None, runFlowerStats=False):
        """Takes the list of flowers for a recursive target and splits them up to fit the given wrapper target(s).
        """
        self.makeChildTargets(flowersAndSizes=runCactusSplitFlowersBySecondaryGrouping(self.flowerNames), 
                              target=target, overlargeTarget=overlargeTarget, phaseNode=phaseNode, runFlowerStats=runFlowerStats)

############################################################
############################################################
############################################################
##The (optional) blast phase that uses the trimming strategy.
############################################################
############################################################
############################################################

def prependUniqueIDs(fas, outputDir):
    """Prepend unique ints to fasta headers.

    (prepend rather than append since trimmed outgroups have a start
    token appended, which complicates removal slightly)
    """
    uniqueID = 0
    ret = []
    for fa in fas:
        outPath = os.path.join(outputDir, os.path.basename(fa))
        out = open(outPath, 'w')
        for line in open(fa):
            if len(line) > 0 and line[0] == '>':
                tokens = line[1:].split()
                tokens[0] = "id=%d|%s" % (uniqueID, tokens[0])
                out.write(">%s\n" % "".join(tokens))
            else:
                out.write(line)
        ret.append(outPath)
        uniqueID += 1
    return ret

def setupDivergenceArgs(cactusWorkflowArguments):
    #Adapt the config file to use arguments for the appropriate divergence distance
    cactusWorkflowArguments.longestPath = getLongestPath(newickTreeParser(cactusWorkflowArguments.speciesTree))
    if cactusWorkflowArguments.outgroupEventNames == None:
        distanceToAddToRootAlignment = getOptionalAttrib(cactusWorkflowArguments.configNode, "distanceToAddToRootAlignment", float, 0.0)
        cactusWorkflowArguments.longestPath += distanceToAddToRootAlignment
    cw = ConfigWrapper(cactusWorkflowArguments.configNode)
    cw.substituteAllDivergenceContolledParametersWithLiterals(cactusWorkflowArguments.longestPath)

def setupFilteringByIdentity(cactusWorkflowArguments):
    #Filter by identity
    cafNode = findRequiredNode(cactusWorkflowArguments.configNode, "caf")
    if getOptionalAttrib(cafNode, "filterByIdentity", bool, False): #Do the identity filtering
        adjustedPath = max(float(cafNode.attrib["identityRatio"]) * cactusWorkflowArguments.longestPath,
        float(cafNode.attrib["minimumDistance"]))
        identity = str(100 - math.ceil(100 * inverseJukesCantor(adjustedPath)))
        cafNode.attrib["lastzArguments"] = cafNode.attrib["lastzArguments"] + (" --identity=%s" % identity)


class CactusTrimmingBlastPhase(CactusPhasesTarget):
    """Blast ingroups vs outgroups using the trimming strategy before
    running cactus setup.
    """
    def run(self):
        # Not worth doing extra work if there aren't any outgroups
        assert self.cactusWorkflowArguments.outgroupEventNames is not None

        self.logToMaster("Running blast using the trimming strategy")

        outgroupsDir = os.path.join(self.getGlobalTempDir(), "outgroupFragments/")
        os.mkdir(outgroupsDir)

        ingroupCoverageDir = os.path.join(self.getGlobalTempDir(), "ingroupCoverageDir/")
        os.mkdir(ingroupCoverageDir)
        self.cactusWorkflowArguments.ingroupCoverageDir = ingroupCoverageDir

        # Get ingroup and outgroup sequences
        exp = ExperimentWrapper(self.cactusWorkflowArguments.experimentNode)
        seqMap = exp.buildSequenceMap()
        # Prepend unique ID to fasta headers to prevent name collision
        renamedInputSeqDir = os.path.join(self.getGlobalTempDir(), "renamedInputs")
        os.mkdir(renamedInputSeqDir)
        uniqueFas = prependUniqueIDs(seqMap.values(), renamedInputSeqDir)
        seqMap = dict(zip(seqMap.keys(), uniqueFas))
        ingroups = map(lambda x: x[1], filter(lambda x: x[0] not in exp.getOutgroupEvents(), seqMap.items()))
        outgroups = [seqMap[i] for i in exp.getOutgroupEvents()]
        self.logToMaster("Ingroup sequences: %s" % (ingroups))
        self.logToMaster("Outgroup sequences: %s" % (outgroups))

        # Change the blast arguments depending on the divergence
        setupDivergenceArgs(self.cactusWorkflowArguments)
        setupFilteringByIdentity(self.cactusWorkflowArguments)

        alignmentsFile = getTempFile("unconvertedAlignments", rootDir=self.getGlobalTempDir())
        findRequiredNode(self.cactusWorkflowArguments.configNode, "caf").attrib["alignments"] = alignmentsFile
        # FIXME: this is really ugly and steals the options from the caf tag
        self.addChildTarget(BlastIngroupsAndOutgroups(
                                          BlastOptions(chunkSize=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "caf"), "chunkSize", int),
                                                        overlapSize=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "caf"), "overlapSize", int),
                                                        lastzArguments=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "caf"), "lastzArguments"),
                                                        compressFiles=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "caf"), "compressFiles", bool),
                                                        realign=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "caf"), "realign", bool), 
                                                        realignArguments=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "caf"), "realignArguments"),
                                                        memory=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "caf"), "lastzMemory", int, sys.maxint),
                                                        minimumSequenceLength=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "caf"), "minimumSequenceLengthForBlast", int, 1),
                                                       trimFlanking=self.getOptionalPhaseAttrib("trimFlanking", int, 10),
                                                       trimMinSize=self.getOptionalPhaseAttrib("trimMinSize", int, 0),
                                                       trimThreshold=self.getOptionalPhaseAttrib("trimThreshold", float, 0.8),
                                                       trimWindowSize=self.getOptionalPhaseAttrib("trimWindowSize", int, 10),
                                                       trimOutgroupFlanking=self.getOptionalPhaseAttrib("trimOutgroupFlanking", int, 100),
                                                       trimOutgroupDepth=self.getOptionalPhaseAttrib("trimOutgroupDepth", int, 1),
                                                       keepParalogs=self.getOptionalPhaseAttrib("keepParalogs", bool, False)), ingroups, outgroups, alignmentsFile, outgroupsDir, ingroupCoverageDir))
        # Point the outgroup sequences to their trimmed versions for
        # phases after this one.
        for outgroup in exp.getOutgroupEvents():
            oldPath = seqMap[outgroup]
            seqMap[outgroup] = os.path.join(outgroupsDir, os.path.basename(oldPath))
        exp.updateTree(exp.getTree(), seqMap)

        self.makeFollowOnPhaseTarget(CactusSetupPhase, "setup")

############################################################
############################################################
############################################################
##The setup phase.
############################################################
############################################################
############################################################

def getLongestPath(node, distance=0.0):

    """Identify the longest path from the mrca of the leaves of the species tree.
    """
    i, j = distance, distance
    if node.left != None:
        i = getLongestPath(node.left, abs(node.left.distance)) + distance
    if node.right != None:  
        j = getLongestPath(node.right, abs(node.right.distance)) + distance
    return max(i, j)

class CactusSetupPhase(CactusPhasesTarget):  
    """Initialises the cactus database and adapts the config file for the run.
    """
    def run(self):
        cw = ConfigWrapper(self.cactusWorkflowArguments.configNode)

        if (not self.cactusWorkflowArguments.configWrapper.getDoTrimStrategy()) or (self.cactusWorkflowArguments.outgroupEventNames == None):
            setupDivergenceArgs(self.cactusWorkflowArguments)

        # we circumvent makeFollowOnPhaseTarget() interface for this job.
        setupTarget = CactusSetupPhase2(cactusWorkflowArguments=self.cactusWorkflowArguments,
                                       phaseName='setup', topFlowerName=self.topFlowerName,
                                       index=0)
        
        #Get the db running and the actual setup going.
        exp = ExperimentWrapper(self.cactusWorkflowArguments.experimentNode)
        if exp.getDbType() == "kyoto_tycoon":
            logger.info("Created ktserver pattern target cactus_setup")
            memory = cw.getKtserverMemory(default=getOptionalAttrib(
                    self.constantsNode, "defaultMemory", int, default=sys.maxint))
            cpu = cw.getKtserverCpu(default=getOptionalAttrib(
                    self.constantsNode, "defaultCpu", int, default=sys.maxint))
            addKtserverDependentChild(self, setupTarget, maxMemory=memory, maxCpu=cpu, isSecondary = False)
        else:
            logger.info("Created follow-on target cactus_setup")
            self.setFollowOnTarget(setupTarget)   
        
class CactusSetupPhase2(CactusPhasesTarget):   
    def run(self):        
        #Now run setup
        messages = runCactusSetup(cactusDiskDatabaseString=self.cactusWorkflowArguments.cactusDiskDatabaseString, 
                       sequences=ExperimentWrapper(self.cactusWorkflowArguments.experimentNode).getSequences(),
                       newickTreeString=self.cactusWorkflowArguments.speciesTree, 
                       outgroupEvents=self.cactusWorkflowArguments.outgroupEventNames,
                       makeEventHeadersAlphaNumeric=self.getOptionalPhaseAttrib("makeEventHeadersAlphaNumeric", bool, False))
        for message in messages:
            self.logToMaster(message)
        self.makeFollowOnPhaseTarget(CactusCafPhase, "caf")
        
############################################################
############################################################
############################################################
#The CAF phase.
#
#Creates the reconstruction structure with blocks
############################################################
############################################################
############################################################

def inverseJukesCantor(d):
    """Takes a substitution distance and calculates the number of expected changes per site (inverse jukes cantor)
    d = -3/4 * log(1 - 4/3 * p)
    exp(-4/3 * d) = 1 - 4/3 * p
    4/3 * p = 1 - exp(-4/3 * d)
    p = 3/4 * (1 - exp(-4/3 * d))
    """
    assert d >= 0.0
    return 0.75 * (1 - math.exp(-d * 4.0/3.0))
    
class CactusCafPhase(CactusPhasesTarget):      
    def run(self):
        if self.cactusWorkflowArguments.ingroupCoverageDir is not None:
            # Convert the bed files to use 64-bit cactus Names instead
            # of the headers. Ideally this should belong in the bar
            # phase but we run stripUniqueIDs before then.
            bedFiles = os.listdir(self.cactusWorkflowArguments.ingroupCoverageDir)
            bedFiles = [os.path.join(self.cactusWorkflowArguments.ingroupCoverageDir, bedFile) for bedFile in bedFiles]
            tempFile = getTempFile(rootDir=self.getGlobalTempDir())
            system("cat %s > %s" % (" ".join(bedFiles), tempFile))
            self.cactusWorkflowArguments.ingroupCoverageFile = getTempFile(rootDir=self.getGlobalTempDir())
            runConvertAlignmentsToInternalNames(self.cactusWorkflowArguments.cactusDiskDatabaseString, tempFile, self.cactusWorkflowArguments.ingroupCoverageFile, self.topFlowerName, isBedFile=True)

        if (not self.cactusWorkflowArguments.configWrapper.getDoTrimStrategy()) or (self.cactusWorkflowArguments.outgroupEventNames == None):
            setupFilteringByIdentity(self.cactusWorkflowArguments)
        #Setup any constraints
        if self.getPhaseIndex() == 0 and self.cactusWorkflowArguments.constraintsFile != None: #Setup the constraints arg
            newConstraintsFile = os.path.join(self.getGlobalTempDir(), "constraints.cig")
            runCactusConvertAlignmentToCactus(self.cactusWorkflowArguments.cactusDiskDatabaseString,
                                              self.cactusWorkflowArguments.constraintsFile, newConstraintsFile)
            self.phaseNode.attrib["constraints"] = newConstraintsFile  
        if self.getOptionalPhaseAttrib("alignments", default="") != "":
            # An alignment file has been provided (likely from the
            # ingroup vs. outgroup blast stage), so just run caf using
            # that file
            assert self.getPhaseNumber() == 1
            convertedAlignmentsFile = getTempFile(rootDir=self.getGlobalTempDir())
            # Convert the cigar file to use 64-bit cactus Names instead of the headers.
            runConvertAlignmentsToInternalNames(self.cactusWorkflowArguments.cactusDiskDatabaseString, self.phaseNode.attrib["alignments"], convertedAlignmentsFile, self.topFlowerName)
            self.logToMaster("Converted headers of cigar file %s to internal names, new file %s" % (self.phaseNode.attrib["alignments"], convertedAlignmentsFile))
            self.phaseNode.attrib["alignments"] = convertedAlignmentsFile
            # While we're at it, remove the unique IDs prepended to
            # the headers inside the cactus DB.
            runStripUniqueIDs(self.cactusWorkflowArguments.cactusDiskDatabaseString)
            self.runPhase(CactusCafWrapperLarge2, CactusBarPhase, "bar")
        elif self.getPhaseIndex()+1 < self.getPhaseNumber(): #Check if there is a repeat phase
            self.runPhase(CactusCafRecursion, CactusCafPhase, "caf", index=self.getPhaseIndex()+1)
        else:
            self.runPhase(CactusCafRecursion, CactusBarPhase, "bar")

class CactusCafRecursion(CactusRecursionTarget):
    """This target does the get flowers down pass for the CAF alignment phase.
    """    
    def run(self):
        self.makeRecursiveTargets()
        self.makeExtendingTargets(target=CactusCafWrapper, overlargeTarget=CactusCafWrapperLarge, runFlowerStats=True)

class CactusCafWrapper(CactusRecursionTarget):
    """Runs cactus_core upon a set of flowers and no alignment file.
    """
    def runCactusCafInWorkflow(self, alignmentFile):
        debugFilePath = self.getOptionalPhaseAttrib("phylogenyDebugPrefix")
        if debugFilePath != None:
            debugFilePath += getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "reference"), "reference")
        messages = runCactusCaf(cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                          alignments=alignmentFile, 
                          flowerNames=self.flowerNames,
                          constraints=self.getOptionalPhaseAttrib("constraints"),  
                          annealingRounds=self.getOptionalPhaseAttrib("annealingRounds"),  
                          deannealingRounds=self.getOptionalPhaseAttrib("deannealingRounds"),
                          trim=self.getOptionalPhaseAttrib("trim"),
                          minimumTreeCoverage=self.getOptionalPhaseAttrib("minimumTreeCoverage", float),
                          blockTrim=self.getOptionalPhaseAttrib("blockTrim", float),
                          minimumBlockDegree=self.getOptionalPhaseAttrib("minimumBlockDegree", int), 
                          minimumIngroupDegree=self.getOptionalPhaseAttrib("minimumIngroupDegree", int),
                          minimumOutgroupDegree=self.getOptionalPhaseAttrib("minimumOutgroupDegree", int),
                          alignmentFilter=self.getOptionalPhaseAttrib("alignmentFilter"),
                          lastzArguments=self.getOptionalPhaseAttrib("lastzArguments"),
                          minimumSequenceLengthForBlast=self.getOptionalPhaseAttrib("minimumSequenceLengthForBlast", int, 1),
                          maxAdjacencyComponentSizeRatio=self.getOptionalPhaseAttrib("maxAdjacencyComponentSizeRatio", float),
                          minLengthForChromosome=self.getOptionalPhaseAttrib("minLengthForChromosome", int),
                          proportionOfUnalignedBasesForNewChromosome=self.getOptionalPhaseAttrib("proportionOfUnalignedBasesForNewChromosome", float),
                          maximumMedianSequenceLengthBetweenLinkedEnds=self.getOptionalPhaseAttrib("maximumMedianSequenceLengthBetweenLinkedEnds", int),
                          realign=self.getOptionalPhaseAttrib("realign", bool),
                          realignArguments=self.getOptionalPhaseAttrib("realignArguments"),
                          phylogenyNumTrees=self.getOptionalPhaseAttrib("phylogenyNumTrees", int, 1),
                          phylogenyRootingMethod=self.getOptionalPhaseAttrib("phylogenyRootingMethod"),
                          phylogenyScoringMethod=self.getOptionalPhaseAttrib("phylogenyScoringMethod"),
                          phylogenyBreakpointScalingFactor=self.getOptionalPhaseAttrib("phylogenyBreakpointScalingFactor"),
                          phylogenySkipSingleCopyBlocks=self.getOptionalPhaseAttrib("phylogenySkipSingleCopyBlocks", bool),
                          phylogenyMaxBaseDistance=self.getOptionalPhaseAttrib("phylogenyMaxBaseDistance"),
                          phylogenyMaxBlockDistance=self.getOptionalPhaseAttrib("phylogenyMaxBlockDistance"),
                          phylogenyDebugFile=debugFilePath,
                          phylogenyKeepSingleDegreeBlocks=self.getOptionalPhaseAttrib("phylogenyKeepSingleDegreeBlocks", bool),
                          phylogenyTreeBuildingMethod=self.getOptionalPhaseAttrib("phylogenyTreeBuildingMethod"),
                          phylogenyCostPerDupPerBase=self.getOptionalPhaseAttrib("phylogenyCostPerDupPerBase"),
                          phylogenyCostPerLossPerBase=self.getOptionalPhaseAttrib("phylogenyCostPerLossPerBase"),
                          referenceEventHeader=getOptionalAttrib(findRequiredNode(self.cactusWorkflowArguments.configNode, "reference"), "reference"),
                          phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce=self.getOptionalPhaseAttrib("phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce"),
                          numTreeBuildingThreads=self.getOptionalPhaseAttrib("numTreeBuildingThreads"),
                          doPhylogeny=self.getOptionalPhaseAttrib("doPhylogeny", bool, True),
                          minimumBlockHomologySupport=self.getOptionalPhaseAttrib("minimumBlockHomologySupport"),
                          minimumBlockDegreeToCheckSupport=self.getOptionalPhaseAttrib("minimumBlockDegreeToCheckSupport"),
                          phylogenyNucleotideScalingFactor=self.getOptionalPhaseAttrib("phylogenyNucleotideScalingFactor"),
                          removeRecoverableChains=self.getOptionalPhaseAttrib("removeRecoverableChains"),
                          minimumNumberOfSpecies=self.getOptionalPhaseAttrib("minimumNumberOfSpecies", int),
                          phylogenyHomologyUnitType=self.getOptionalPhaseAttrib("phylogenyHomologyUnitType"),
                          phylogenyDistanceCorrectionMethod=self.getOptionalPhaseAttrib("phylogenyDistanceCorrectionMethod"),
                          maxRecoverableChainsIterations=self.getOptionalPhaseAttrib("maxRecoverableChainsIterations", int),
                          maxRecoverableChainLength=self.getOptionalPhaseAttrib("maxRecoverableChainLength", int))
        for message in messages:
            self.logToMaster(message)
    
    def run(self):
        self.runCactusCafInWorkflow(alignmentFile=None)
       
class CactusCafWrapperLarge(CactusRecursionTarget):
    """Runs blast on the given flower and passes the resulting alignment to cactus core.
    """
    def run(self):
        logger.info("Starting the cactus aligner target")
        #Generate a temporary file to hold the alignments
        alignmentFile = os.path.join(self.getGlobalTempDir(), "alignments.cigar")
        flowerName = decodeFirstFlowerName(self.flowerNames)
        self.addChildTarget(BlastFlower(self.cactusDiskDatabaseString, 
                                          flowerName, alignmentFile, 
                                          blastOptions=\
                                          BlastOptions(chunkSize=self.getOptionalPhaseAttrib("chunkSize", int),
                                                        overlapSize=self.getOptionalPhaseAttrib("overlapSize", int),
                                                        lastzArguments=self.getOptionalPhaseAttrib("lastzArguments"),
                                                        compressFiles=self.getOptionalPhaseAttrib("compressFiles", bool),
                                                        realign=self.getOptionalPhaseAttrib("realign", bool), 
                                                        realignArguments=self.getOptionalPhaseAttrib("realignArguments"),
                                                        memory=self.getOptionalPhaseAttrib("lastzMemory", int, sys.maxint),
                                                        minimumSequenceLength=self.getOptionalPhaseAttrib("minimumSequenceLengthForBlast", int, 1))))
        #Now setup a call to cactus core wrapper as a follow on
        self.phaseNode.attrib["alignments"] = alignmentFile
        self.makeFollowOnRecursiveTarget(CactusCafWrapperLarge2)
        
class CactusCafWrapperLarge2(CactusCafWrapper):
    """Runs cactus_core upon a one flower and one alignment file.
    """
    def run(self):
        self.runCactusCafInWorkflow(alignmentFile=self.phaseNode.attrib["alignments"])
        
############################################################
############################################################
############################################################
#The BAR phase.
#
#Creates the reconstruction structure with blocks
############################################################
############################################################
############################################################

class CactusBarPhase(CactusPhasesTarget): 
    """Runs bar algorithm
    """  
    def run(self):
        self.runPhase(CactusBarRecursion, CactusNormalPhase, "normal", doRecursion=self.getOptionalPhaseAttrib("runBar", bool, False))

class CactusBarRecursion(CactusRecursionTarget):
    """This target does the get flowers down pass for the BAR alignment phase.
    """
    def run(self):
        self.makeRecursiveTargets()
        self.makeExtendingTargets(target=CactusBarWrapper, overlargeTarget=CactusBarWrapperLarge, runFlowerStats=True)

def runBarForTarget(self, calculateWhichEndsToComputeSeparately=None, endAlignmentsToPrecomputeOutputFile=None, precomputedAlignments=None):
    return runCactusBar(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                 flowerNames=self.flowerNames, 
                 maximumLength=self.getOptionalPhaseAttrib("bandingLimit", float),
                 spanningTrees=self.getOptionalPhaseAttrib("spanningTrees", int), 
                 gapGamma=self.getOptionalPhaseAttrib( "gapGamma", float), 
                 matchGamma=self.getOptionalPhaseAttrib( "matchGamma", float), 
                 splitMatrixBiggerThanThis=self.getOptionalPhaseAttrib("splitMatrixBiggerThanThis", int), 
                 anchorMatrixBiggerThanThis=self.getOptionalPhaseAttrib("anchorMatrixBiggerThanThis", int), 
                 repeatMaskMatrixBiggerThanThis=self.getOptionalPhaseAttrib("repeatMaskMatrixBiggerThanThis", int), 
                 diagonalExpansion=self.getOptionalPhaseAttrib("diagonalExpansion"),
                 constraintDiagonalTrim=self.getOptionalPhaseAttrib("constraintDiagonalTrim", int), 
                 minimumBlockDegree=self.getOptionalPhaseAttrib("minimumBlockDegree", int),
                 minimumIngroupDegree=self.getOptionalPhaseAttrib("minimumIngroupDegree", int),
                 minimumOutgroupDegree=self.getOptionalPhaseAttrib("minimumOutgroupDegree", int),
                 alignAmbiguityCharacters=self.getOptionalPhaseAttrib("alignAmbiguityCharacters", bool),
                 pruneOutStubAlignments=self.getOptionalPhaseAttrib("pruneOutStubAlignments", bool),
                 useProgressiveMerging=self.getOptionalPhaseAttrib("useProgressiveMerging", bool),
                 calculateWhichEndsToComputeSeparately=calculateWhichEndsToComputeSeparately,
                 endAlignmentsToPrecomputeOutputFile=endAlignmentsToPrecomputeOutputFile,
                 largeEndSize=self.getOptionalPhaseAttrib("largeEndSize", int),
                 precomputedAlignments=precomputedAlignments,
                 ingroupCoverageFile=self.cactusWorkflowArguments.ingroupCoverageFile if self.getOptionalPhaseAttrib("rescue", bool) else None,
                 minimumSizeToRescue=self.getOptionalPhaseAttrib("minimumSizeToRescue"),
                 minimumCoverageToRescue=self.getOptionalPhaseAttrib("minimumCoverageToRescue"),
                 minimumNumberOfSpecies=self.getOptionalPhaseAttrib("minimumNumberOfSpecies", int))

class CactusBarWrapper(CactusRecursionTarget):
    """Runs the BAR algorithm implementation.
    """
    def run(self):
        messages = runBarForTarget(self)
        for message in messages:
            self.logToMaster(message)       
        
class CactusBarWrapperLarge(CactusRecursionTarget):
    """Breaks up the bar into a series of smaller bars, then runs them.
    """
    def run(self):
        logger.info("Starting the cactus bar preprocessor target to breakup the bar alignment")
        precomputedAlignmentFiles = []
        veryLargeEndSize=self.getOptionalPhaseAttrib("veryLargeEndSize", int, default=1000000)
        maxFlowerGroupSize = self.getOptionalTargetAttrib("maxFlowerGroupSize", int, 
                                            default=CactusRecursionTarget.maxSequenceSizeOfFlowerGroupingDefault)
        endsToAlign = []
        totalSize = 0
        alignmentFileCount = 0
        for line in runBarForTarget(self, calculateWhichEndsToComputeSeparately=True):
            endToAlign, sequencesInEndAlignment, basesInEndAlignment = line.split()
            sequencesInEndAlignment = int(sequencesInEndAlignment)
            basesInEndAlignment = int(basesInEndAlignment)
            #If we have a really big end align separately
            if basesInEndAlignment >= veryLargeEndSize:
                self.addChildTarget(CactusBarEndAlignerWrapper(self.phaseNode, self.constantsNode, self.cactusDiskDatabaseString, self.flowerNames, 
                                                           True, self.cactusWorkflowArguments, [ endToAlign ], os.path.join(self.getGlobalTempDir(), "endAlignments.%i" % alignmentFileCount)))
                self.logToMaster("Precomputing very large end alignment for %s with %i caps and %i bases" % \
                             (endToAlign, sequencesInEndAlignment, basesInEndAlignment))
                alignmentFileCount += 1
            else:
                endsToAlign.append(endToAlign)
                totalSize += basesInEndAlignment
                if totalSize >= maxFlowerGroupSize:
                    self.addChildTarget(CactusBarEndAlignerWrapper(self.phaseNode, self.constantsNode, self.cactusDiskDatabaseString, self.flowerNames, 
                                                           False, self.cactusWorkflowArguments, endsToAlign, os.path.join(self.getGlobalTempDir(), "endAlignments.%i" % alignmentFileCount)))
                    endsToAlign = []
                    totalSize = 0
                    alignmentFileCount += 1
        if len(endsToAlign) > 0:
            self.addChildTarget(CactusBarEndAlignerWrapper(self.phaseNode, self.constantsNode, self.cactusDiskDatabaseString, self.flowerNames, 
                                                           False, self.cactusWorkflowArguments, endsToAlign, os.path.join(self.getGlobalTempDir(), "endAlignments.%i" % alignmentFileCount)))
            alignmentFileCount += 1
        self.phaseNode.attrib["precomputedAlignmentFiles"] = " ".join([ os.path.join(self.getGlobalTempDir(), ("endAlignments.%i") % i) for i in range(alignmentFileCount) ]) 
        self.makeFollowOnRecursiveTarget(CactusBarWrapperWithPrecomputedEndAlignments)
        self.logToMaster("Breaking bar job into %i separate jobs" % \
                             (alignmentFileCount))

class CactusBarEndAlignerWrapper(CactusRecursionTarget):
    """Computes an end alignment.
    """
    def __init__(self, phaseNode, constantsNode, cactusDiskDatabaseString, flowerNames, overlarge, cactusWorkflowArguments, endsToAlign, alignmentFile):
        CactusRecursionTarget.__init__(self, phaseNode, constantsNode, cactusDiskDatabaseString, flowerNames, overlarge, cactusWorkflowArguments)
        self.endsToAlign = endsToAlign
        self.alignmentFile = alignmentFile
    
    def run(self):
        self.endsToAlign = [ int(i) for i in self.endsToAlign ]
        self.endsToAlign.sort()
        self.flowerNames = encodeFlowerNames((decodeFirstFlowerName(self.flowerNames),) + tuple(self.endsToAlign)) #The ends to align become like extra flower names
        messages = runBarForTarget(self, 
                                   endAlignmentsToPrecomputeOutputFile=self.alignmentFile)
        for message in messages:
            self.logToMaster(message)
        
class CactusBarWrapperWithPrecomputedEndAlignments(CactusRecursionTarget):
    """Runs the BAR algorithm implementation with some precomputed end alignments.
    """
    def run(self):
        if self.phaseNode.attrib["precomputedAlignmentFiles"] != "":
            messages = runBarForTarget(self, precomputedAlignments=self.phaseNode.attrib["precomputedAlignmentFiles"])
        else:
            messages = runBarForTarget(self)
        for message in messages:
            self.logToMaster(message)
        
############################################################
############################################################
############################################################
#Normalisation pass
############################################################
############################################################
############################################################
    
class CactusNormalPhase(CactusPhasesTarget):
    """Phase to normalise the graph, ensuring all chains are maximal
    """
    def run(self):
        normalisationIterations = self.getOptionalPhaseAttrib("iterations", int, default=0)
        if normalisationIterations > 0:
            self.phaseNode.attrib["normalised"] = "1"
            self.phaseNode.attrib["iterations"] = str(normalisationIterations-1)
            self.runPhase(CactusNormalRecursion, CactusNormalPhase, "normal")
        else:
            self.makeFollowOnPhaseTarget(CactusAVGPhase, "avg")
     
class CactusNormalRecursion(CactusRecursionTarget):
    """This target does the down pass for the normal phase.
    """
    def run(self):
        self.makeRecursiveTargets()
        self.makeFollowOnRecursiveTarget(CactusNormalRecursion2)
        
class CactusNormalRecursion2(CactusRecursionTarget):
    """This target sets up the normal wrapper in an up traversal of the tree.
    """
    def run(self):
        self.makeWrapperTargets(CactusNormalWrapper)
        
class CactusNormalWrapper(CactusRecursionTarget):
    """This targets run the normalisation script.
    """ 
    def run(self):
        runCactusMakeNormal(self.cactusDiskDatabaseString, flowerNames=self.flowerNames, 
                            maxNumberOfChains=self.getOptionalPhaseAttrib("maxNumberOfChains", int, default=30))

############################################################
############################################################
############################################################
#Phylogeny pass
############################################################
############################################################
############################################################
    
class CactusAVGPhase(CactusPhasesTarget): 
    """Phase to build avgs for each flower.
    """       
    def run(self):
        self.runPhase(CactusAVGRecursion, CactusReferencePhase, "reference", doRecursion=self.getOptionalPhaseAttrib("buildAvgs", bool, False))

class CactusAVGRecursion(CactusRecursionTarget):
    """This target does the recursive pass for the AVG phase.
    """
    def run(self):
        self.makeFollowOnRecursiveTarget(CactusAVGRecursion2)
        self.makeWrapperTargets(CactusAVGWrapper)

class CactusAVGRecursion2(CactusRecursionTarget):
    """This target does the recursive pass for the AVG phase.
    """
    def run(self):
        self.makeRecursiveTargets(target=CactusAVGRecursion)

class CactusAVGWrapper(CactusRecursionTarget):
    """This target runs tree building
    """
    def run(self):
        runCactusPhylogeny(self.cactusDiskDatabaseString, flowerNames=self.flowerNames)

############################################################
############################################################
############################################################
#Reference pass
############################################################
############################################################
############################################################
    
class CactusReferencePhase(CactusPhasesTarget):     
    def run(self):
        """Runs the reference problem algorithm
        """
        self.setupSecondaryDatabase()
        self.phaseNode.attrib["experimentPath"] = self.cactusWorkflowArguments.experimentFile
        self.phaseNode.attrib["secondaryDatabaseString"] = self.cactusWorkflowArguments.secondaryDatabaseString
        self.runPhase(CactusReferenceRecursion, CactusSetReferenceCoordinatesDownPhase, "reference", 
                      doRecursion=self.getOptionalPhaseAttrib("buildReference", bool, False),
                      launchSecondaryKtForRecursiveTarget = True)
        
class CactusReferenceRecursion(CactusRecursionTarget):
    """This target creates the wrappers to run the reference problem algorithm, the follow on target then recurses down.
    """
    def run(self):
        self.makeWrapperTargets(CactusReferenceWrapper, runFlowerStats=True)
        self.makeFollowOnRecursiveTarget(CactusReferenceRecursion2)
        
class CactusReferenceWrapper(CactusRecursionTarget):
    """Actually run the reference code.
    """
    def run(self):
        for message in runCactusReference(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                       flowerNames=self.flowerNames, 
                       matchingAlgorithm=self.getOptionalPhaseAttrib("matchingAlgorithm"), 
                       permutations=self.getOptionalPhaseAttrib("permutations", int),
                       referenceEventString=self.getOptionalPhaseAttrib("reference"),
                       useSimulatedAnnealing=self.getOptionalPhaseAttrib("useSimulatedAnnealing", bool),
                       theta=self.getOptionalPhaseAttrib("theta", float),
                       phi=self.getOptionalPhaseAttrib("phi", float),
                       maxWalkForCalculatingZ=self.getOptionalPhaseAttrib("maxWalkForCalculatingZ", int),
                       ignoreUnalignedGaps=self.getOptionalPhaseAttrib("ignoreUnalignedGaps", bool),
                       wiggle=self.getOptionalPhaseAttrib("wiggle", float),
                       numberOfNs=self.getOptionalPhaseAttrib("numberOfNs", int),
                       minNumberOfSequencesToSupportAdjacency=self.getOptionalPhaseAttrib("minNumberOfSequencesToSupportAdjacency", int),
                       makeScaffolds=self.getOptionalPhaseAttrib("makeScaffolds", bool)):
            self.logToMaster(message)

class CactusReferenceRecursion2(CactusRecursionTarget):
    def run(self):
        self.makeRecursiveTargets(target=CactusReferenceRecursion)
        self.makeFollowOnRecursiveTarget(CactusReferenceRecursion3)
        
class CactusReferenceRecursion3(CactusRecursionTarget):
    """After completing the recursion for the reference algorithm, the up pass of adding in the reference coordinates is performed.
    """
    def run(self):
        self.makeWrapperTargets(CactusSetReferenceCoordinatesUpWrapper)

class CactusSetReferenceCoordinatesUpWrapper(CactusRecursionTarget):
    """Does the up pass for filling in the reference sequence coordinates, once a reference has been established.
    """ 
    def run(self):
        runCactusAddReferenceCoordinates(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                                         secondaryDatabaseString=self.getOptionalPhaseAttrib("secondaryDatabaseString"),
                                         flowerNames=self.flowerNames,
                                         referenceEventString=self.getOptionalPhaseAttrib("reference"), 
                                         outgroupEventString=self.getOptionalPhaseAttrib("outgroup"), 
                                         bottomUpPhase=True)
        
class CactusSetReferenceCoordinatesDownPhase(CactusPhasesTarget):
    """This is the second part of the reference coordinate setting, the down pass.
    """
    def run(self):
        self.cleanupSecondaryDatabase()
        self.runPhase(CactusSetReferenceCoordinatesDownRecursion, CactusExtractReferencePhase, "check", doRecursion=self.getOptionalPhaseAttrib("buildReference", bool, False))
        
class CactusSetReferenceCoordinatesDownRecursion(CactusRecursionTarget):
    """Does the down pass for filling Fills in the coordinates, once a reference is added.
    """        
    def run(self):
        self.makeWrapperTargets(CactusSetReferenceCoordinatesDownWrapper)
        self.makeFollowOnRecursiveTarget(CactusSetReferenceCoordinatesDownRecursion2)

class CactusSetReferenceCoordinatesDownRecursion2(CactusRecursionTarget):
    def run(self):
        self.makeRecursiveTargets(target=CactusSetReferenceCoordinatesDownRecursion)
        
class CactusSetReferenceCoordinatesDownWrapper(CactusRecursionTarget):
    """Does the down pass for filling Fills in the coordinates, once a reference is added.
    """        
    def run(self):
        runCactusAddReferenceCoordinates(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                                         flowerNames=self.flowerNames,
                                         referenceEventString=self.getOptionalPhaseAttrib("reference"),
                                         outgroupEventString=self.getOptionalPhaseAttrib("outgroup"), 
                                         bottomUpPhase=False)

class CactusExtractReferencePhase(CactusPhasesTarget):
    def run(self):
        if hasattr(self.cactusWorkflowArguments, 'buildReference') and\
               self.cactusWorkflowArguments.buildReference:
            self.logToMaster("Starting Reference Extract Phase")
            experiment = ExperimentWrapper(self.cactusWorkflowArguments.experimentNode)
            if experiment.getReferencePath() is not None:
                eventName = os.path.basename(experiment.getReferencePath())
                if eventName.find('.') >= 0:
                    eventName = eventName[:eventName.rfind('.')]
                    cmdLine = "cactus_getReferenceSeq --cactusDisk \'%s\' --flowerName 0 --referenceEventString %s --outputFile %s --logLevel %s" % \
                              (experiment.getDiskDatabaseString(), eventName,
                               experiment.getReferencePath(), getLogLevelString())                        
                    system(cmdLine)          
        self.makeFollowOnPhaseTarget(CactusCheckPhase, "check")

############################################################
############################################################
############################################################
#Check pass
############################################################
############################################################
############################################################
    
class CactusCheckPhase(CactusPhasesTarget):
    """The check phase, where we verify everything is as it should be
    """
    def run(self):
        normalNode = findRequiredNode(self.cactusWorkflowArguments.configNode, "normal")
        self.phaseNode.attrib["checkNormalised"] = getOptionalAttrib(normalNode, "normalised", default="0")
        self.runPhase(CactusCheckRecursion, CactusHalGeneratorPhase, "hal", doRecursion=self.getOptionalPhaseAttrib("runCheck", bool, False))
        
class CactusCheckRecursion(CactusRecursionTarget):
    """This target does the recursive pass for the check phase.
    """
    def run(self):
        self.makeRecursiveTargets()
        self.makeWrapperTargets(CactusCheckWrapper)
        
class CactusCheckWrapper(CactusRecursionTarget):
    """Runs the actual check wrapper
    """
    def run(self):
        runCactusCheck(self.cactusDiskDatabaseString, self.flowerNames, checkNormalised=self.getOptionalPhaseAttrib("checkNormalised", bool, False))

############################################################
############################################################
############################################################
#Hal generation
############################################################
############################################################
############################################################

class CactusHalGeneratorPhase(CactusPhasesTarget):
    def run(self):
        referenceNode = findRequiredNode(self.cactusWorkflowArguments.configNode, "reference")
        if referenceNode.attrib.has_key("reference"):
            self.phaseNode.attrib["reference"] = referenceNode.attrib["reference"]
        if self.getOptionalPhaseAttrib("buildFasta", bool, default=False):
            self.phaseNode.attrib["fastaPath"] = self.cactusWorkflowArguments.experimentNode.find("hal").attrib["fastaPath"]
            self.makeRecursiveChildTarget(CactusFastaGenerator)
        if self.getOptionalPhaseAttrib("buildHal", bool, default=False):
            self.setupSecondaryDatabase()
            self.phaseNode.attrib["experimentPath"] = self.cactusWorkflowArguments.experimentFile
            self.phaseNode.attrib["secondaryDatabaseString"] = self.cactusWorkflowArguments.secondaryDatabaseString
            self.phaseNode.attrib["outputFile"]=self.cactusWorkflowArguments.experimentNode.find("hal").attrib["halPath"]
            self.makeFollowOnPhaseTarget(CactusHalGeneratorPhase2, "hal")
            self.makeRecursiveChildTarget(CactusHalGeneratorRecursion, launchSecondaryKtForRecursiveTarget=True)

class CactusFastaGenerator(CactusRecursionTarget):
    def run(self):
        runCactusFastaGenerator(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                                    flowerName=decodeFirstFlowerName(self.flowerNames),
                                    outputFile=self.getOptionalPhaseAttrib("fastaPath"),
                                    referenceEventString=self.getOptionalPhaseAttrib("reference"))
            
class CactusHalGeneratorPhase2(CactusHalGeneratorPhase):
    def run(self): 
        self.cleanupSecondaryDatabase()

class CactusHalGeneratorRecursion(CactusRecursionTarget):
    """Generate the hal file by merging indexed hal files from the children.
    """ 
    def run(self):
        i = extractNode(self.phaseNode)
        if "outputFile" in i.attrib:
            i.attrib.pop("outputFile")
        self.makeRecursiveTargets(phaseNode=i)
        self.makeFollowOnRecursiveTarget(CactusHalGeneratorUpWrapper)

class CactusHalGeneratorUpWrapper(CactusRecursionTarget):
    """Does the up pass for filling in the coordinates, once a reference is added.
    """ 
    def run(self):
        runCactusHalGenerator(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                              secondaryDatabaseString=self.getOptionalPhaseAttrib("secondaryDatabaseString"),
                              flowerNames=self.flowerNames,
                              referenceEventString=self.getOptionalPhaseAttrib("reference"), #self.configNode.attrib["reference"], #self.getOptionalPhaseAttrib("reference"), 
                              outputFile=self.getOptionalPhaseAttrib("outputFile"),
                              showOnlySubstitutionsWithRespectToReference=\
                              self.getOptionalPhaseAttrib("showOnlySubstitutionsWithRespectToReference", bool))

class CactusHalGeneratorPhaseCleanup(CactusPhasesTarget):
    """Cleanup the database used to build the hal
    """
    def run(self):
        self.cleanupSecondaryDatabase()

############################################################
############################################################
############################################################
#Main function
############################################################
############################################################
############################################################
class CactusWorkflowArguments:
    """Object for representing a cactus workflow's arguments
    """
    def __init__(self, options):
        self.experimentFile = getTempFile("tempExperimentFileCopy", rootDir=os.path.dirname(options.experimentFile))
        shutil.copyfile(options.experimentFile, self.experimentFile)
        self.experimentNode = ET.parse(self.experimentFile).getroot()
        self.experimentWrapper = ExperimentWrapper(self.experimentNode)
        #Get the database string
        self.cactusDiskDatabaseString = ET.tostring(self.experimentNode.find("cactus_disk").find("st_kv_database_conf")).translate(None, '\n')
        #Get the species tree
        self.speciesTree = self.experimentNode.attrib["species_tree"]
        #Get any list of 'required species' for the blocks of the cactus.
        self.outgroupEventNames = getOptionalAttrib(self.experimentNode, "outgroup_events")
        #Constraints
        self.constraintsFile = getOptionalAttrib(self.experimentNode, "constraints")
        #Space to put the path to the directory containing beds of
        #outgroup coverage on ingroups, so that any sequence aligning
        #to an outgroup can be rescued after bar phase
        self.ingroupCoverageDir = None
        # Same, but for the final bed file
        self.ingroupCoverageFile = None
        #Secondary, scratch DB
        secondaryConf = copy.deepcopy(self.experimentNode.find("cactus_disk").find("st_kv_database_conf"))
        secondaryElem = DbElemWrapper(secondaryConf)
        dbPath = secondaryElem.getDbDir()
        assert dbPath is not None
        secondaryDbPath = os.path.join(os.path.dirname(dbPath), "%s_tempSecondaryDatabaseDir_%s" % (
            os.path.basename(dbPath), random.random()))
        secondaryElem.setDbDir(secondaryDbPath)
        if secondaryElem.getDbType() == "kyoto_tycoon":
            secondaryElem.setDbPort(secondaryElem.getDbPort() + 100)
        self.secondaryDatabaseString = secondaryElem.getConfString()
            
        #The config node
        self.configNode = ET.parse(self.experimentWrapper.getConfigPath()).getroot()
        self.configWrapper = ConfigWrapper(self.configNode)
        #Now deal with the constants that ned to be added here
        self.configWrapper.substituteAllPredefinedConstantsWithLiterals()
        self.configWrapper.setBuildHal(options.buildHal)
        self.configWrapper.setBuildFasta(options.buildFasta)
        
        #Now build the remaining options from the arguments
        if options.buildAvgs:
            findRequiredNode(self.configNode, "avg").attrib["buildAvgs"] = "1"
        if options.buildReference:
            findRequiredNode(self.configNode, "reference").attrib["buildReference"] = "1"
            

def addCactusWorkflowOptions(parser):
    parser.add_option("--experiment", dest="experimentFile", 
                      help="The file containing a link to the experiment parameters")
    
    parser.add_option("--buildAvgs", dest="buildAvgs", action="store_true",
                      help="Build trees", default=False)
    
    parser.add_option("--buildReference", dest="buildReference", action="store_true",
                      help="Creates a reference ordering for the flowers", default=False)
    
    parser.add_option("--buildHal", dest="buildHal", action="store_true",
                      help="Build a hal file", default=False)
    
    parser.add_option("--buildFasta", dest="buildFasta", action="store_true",
                      help="Build a fasta file of the input sequences (and reference sequence, used with hal output)", 
                      default=False)

    parser.add_option("--test", dest="test", action="store_true",
                      help="Run doctest unit tests")

class RunCactusPreprocessorThenCactusSetup(Target):
    def __init__(self, options):
        Target.__init__(self)
        self.options = options
        
    def run(self):
        cactusWorkflowArguments=CactusWorkflowArguments(self.options)
        eW = ExperimentWrapper(cactusWorkflowArguments.experimentNode)
        outputSequenceFiles = CactusPreprocessor.getOutputSequenceFiles(eW.getSequences(), eW.getOutputSequenceDir())
        self.addChildTarget(CactusPreprocessor(eW.getSequences(), outputSequenceFiles, cactusWorkflowArguments.configNode))
        #Now make the setup, replacing the input sequences with the preprocessed sequences
        eW.setSequences(outputSequenceFiles)
        self.logToMaster("doTrimStrategy() = %s, outgroupEventNames = %s" % (cactusWorkflowArguments.configWrapper.getDoTrimStrategy(), cactusWorkflowArguments.outgroupEventNames))
        if cactusWorkflowArguments.configWrapper.getDoTrimStrategy() and cactusWorkflowArguments.outgroupEventNames is not None:
            # Use the trimming strategy to blast ingroups vs outgroups.
            self.setFollowOnTarget(CactusTrimmingBlastPhase(cactusWorkflowArguments=cactusWorkflowArguments, phaseName="trimBlast"))
        else:
            self.setFollowOnTarget(CactusSetupPhase(cactusWorkflowArguments=cactusWorkflowArguments,
                                                    phaseName="setup"))
        
def main():
    ##########################################
    #Construct the arguments.
    ##########################################
    
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)
    addCactusWorkflowOptions(parser)
        
    options, args = parser.parse_args()
    if options.test:
        _test()
    setLoggingFromOptions(options)
    
    if len(args) != 0:
        raise RuntimeError("Unrecognised input arguments: %s" % " ".join(args))

    cactusWorkflowArguments = CactusWorkflowArguments(options)
    Stack(RunCactusPreprocessorThenCactusSetup(options)).startJobTree(options)

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.pipeline.cactus_workflow import *
    main()
