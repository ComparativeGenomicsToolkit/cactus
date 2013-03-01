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
from optparse import OptionParser

from sonLib.bioio import getTempFile
from sonLib.bioio import newickTreeParser

from sonLib.bioio import logger
from sonLib.bioio import setLoggingFromOptions
from sonLib.bioio import getTempDirectory
from sonLib.bioio import system

from cactus.shared.common import cactusRootPath
  
from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import getLogLevelString

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

from cactus.blast.cactus_blast import BlastFlower
from cactus.blast.cactus_blast import BlastOptions

from cactus.preprocessor.cactus_preprocessor import BatchPreprocessor
from cactus.preprocessor.cactus_preprocessor import PreprocessorHelper

from cactus.progressive.experimentWrapper import ExperimentWrapper
from cactus.progressive.experimentWrapper import DbElemWrapper
from cactus.pipeline.ktserverJobTree import addKtserverDependentChild

############################################################
############################################################
############################################################
##Shared functions
############################################################
############################################################
############################################################

def getOptionalAttrib(node, attribName, typeFn=None, default=None):
    """Get an optional attrib, or None, if not set or node is None
    """
    if node != None and node.attrib.has_key(attribName):
        if typeFn != None:
            if typeFn == bool:
                return bool(int(node.attrib[attribName]))
            return typeFn(node.attrib[attribName])
        return node.attrib[attribName]
    return default

def findRequiredNode(configNode, nodeName, index=0):
    """Retrieve an xml node, complain if its not there.
    """
    nodes = configNode.findall(nodeName)
    if nodes == None:
        raise RuntimeError("Could not find any nodes with name %s in %s node" % (nodeName, configNode))
    if index >= len(nodes):
        raise RuntimeError("Could not find a node with name %s and index %i in %s node" % (nodeName, index, configNode))
    return nodes[index]

def extractNode(node):
    """Make an XML node free of its parent subtree
    """
    return ET.fromstring(ET.tostring(node))

def getTargetNode(phaseNode, targetClass):
    """Gets a target node for a given target.
    """
    className = str(targetClass).split(".")[-1]
    assert className != ''
    return phaseNode.find(className)

class CactusTarget(Target):
    """Base target for all cactus workflow targets.
    """
    def __init__(self, phaseNode, overlarge=False):
        self.phaseNode = phaseNode
        self.overlarge = overlarge
        self.targetNode = getTargetNode(self.phaseNode, self.__class__)
        if overlarge:
            Target.__init__(self, memory=self.getOptionalTargetAttrib("overlargeMemory", typeFn=int, default=sys.maxint), 
                            cpu=self.getOptionalTargetAttrib("overlargeCpu", typeFn=int, default=sys.maxint))
        else:
            Target.__init__(self, memory=self.getOptionalTargetAttrib("memory", typeFn=int, default=sys.maxint), 
                            cpu=self.getOptionalTargetAttrib("cpu", typeFn=int, default=sys.maxint))
    
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
        CactusTarget.__init__(self, phaseNode=phaseNode, overlarge=False)
        self.index = index
        self.cactusWorkflowArguments = cactusWorkflowArguments
        self.topFlowerName = topFlowerName
    
    def makeRecursiveChildTarget(self, target, launchSecondaryKtForRecursiveTarget=False):
        newChild = target(phaseNode=extractNode(self.phaseNode), 
                          cactusDiskDatabaseString=self.cactusWorkflowArguments.cactusDiskDatabaseString, 
                          flowerNames=encodeFlowerNames((self.topFlowerName,)), overlarge=True)
        
        if launchSecondaryKtForRecursiveTarget and ExperimentWrapper(self.cactusWorkflowArguments.experimentNode).getDbType() == "kyoto_tycoon":
            addKtserverDependentChild(self, newChild, isSecondary = True)
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
        # if kyoto_tycoon is true, then the database will be created by the call to
        # addKtserverDependentChild() within makeRecursiveChildTarget()
    
    def cleanupSecondaryDatabase(self):
        """Cleanup the secondary database
        """
        confXML = ET.fromstring(self.cactusWorkflowArguments.secondaryDatabaseString)
        dbElem = DbElemWrapper(confXML)
        if dbElem.getDbType() != "kyoto_tycoon":
            runCactusSecondaryDatabase(self.cactusWorkflowArguments.secondaryDatabaseString, create=False)
        # if kyoto_tycoon is true, then the database will be killed by the call to
        # addKtserverDependentChild() within makeRecursiveChildTarget()

class CactusRecursionTarget(CactusTarget):
    """Base recursive target for traversals up and down the cactus tree.
    """
    maxSequenceSizeOfFlowerGroupingDefault = 1000000
    def __init__(self, phaseNode, cactusDiskDatabaseString, flowerNames, overlarge=False):
        CactusTarget.__init__(self, phaseNode=phaseNode, overlarge=overlarge)
        self.cactusDiskDatabaseString = cactusDiskDatabaseString
        self.flowerNames = flowerNames  
        
    def makeFollowOnRecursiveTarget(self, target, phaseNode=None):
        """Sets the followon to the given recursive target
        """
        if phaseNode == None:
            phaseNode = self.phaseNode
        self.setFollowOnTarget(target(phaseNode=phaseNode, 
                                   cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                                   flowerNames=self.flowerNames, overlarge=self.overlarge))
        
    def makeChildTargets(self, flowersAndSizes, target, overlargeTarget=None, 
                         phaseNode=None):
        """Make a set of child targets for a given set of flowers and chosen child target
        """
        if overlargeTarget == None:
            overlargeTarget = target
        if phaseNode == None:
            phaseNode = self.phaseNode
        for overlarge, flowerNames in flowersAndSizes:
            if overlarge: #Make sure large flowers are on there own, in their own job
                flowerStatsString = runCactusFlowerStats(cactusDiskDatabaseString=self.cactusDiskDatabaseString, flowerName=decodeFirstFlowerName(flowerNames))
                self.logToMaster("Adding an oversize flower for target class %s and stats %s" \
                                         % (overlargeTarget, flowerStatsString))
                self.addChildTarget(overlargeTarget(cactusDiskDatabaseString=self.cactusDiskDatabaseString, phaseNode=phaseNode, 
                                                    flowerNames=flowerNames, overlarge=True)) #This ensures overlarge flowers, 
            else:
                self.addChildTarget(target(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                                           phaseNode=phaseNode, flowerNames=flowerNames, overlarge=False))
        
    def makeRecursiveTargets(self, target=None, phaseNode=None):
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
                              target=target, phaseNode=phaseNode)
    
    def makeExtendingTargets(self, target, overlargeTarget=None, phaseNode=None):
        """Make set of child targets that extend the current cactus tree.
        """
        targetNode = getTargetNode(self.phaseNode, target)
        flowersAndSizes=runCactusExtendFlowers(cactusDiskDatabaseString=self.cactusDiskDatabaseString, flowerNames=self.flowerNames, 
                                              minSequenceSizeOfFlower=getOptionalAttrib(targetNode, "minFlowerSize", int, 1), 
                                              maxSequenceSizeOfFlowerGrouping=getOptionalAttrib(targetNode, "maxFlowerGroupSize", int, 
                                              default=CactusRecursionTarget.maxSequenceSizeOfFlowerGroupingDefault))
        self.makeChildTargets(flowersAndSizes=flowersAndSizes, 
                              target=target, overlargeTarget=overlargeTarget, 
                              phaseNode=phaseNode)
    
    def makeWrapperTargets(self, target, overlargeTarget=None, phaseNode=None):
        """Takes the list of flowers for a recursive target and splits them up to fit the given wrapper target(s).
        """
        self.makeChildTargets(flowersAndSizes=runCactusSplitFlowersBySecondaryGrouping(self.flowerNames), 
                              target=target, overlargeTarget=overlargeTarget, phaseNode=phaseNode)
        
############################################################
############################################################
############################################################
##The preprocessor phase.
############################################################
############################################################
############################################################

class CactusPreprocessorPhase(CactusPhasesTarget):
    def run(self):
        self.logToMaster("Starting preprocessor phase target at %s seconds" % time.time())
        tempDir = getTempDirectory(self.getGlobalTempDir())
        prepHelper = PreprocessorHelper(self.cactusWorkflowArguments, self.cactusWorkflowArguments.sequences)
        processedSequences = []
        for sequence in self.cactusWorkflowArguments.sequences:
            prepXmlElems = prepHelper.getFilteredXmlElems(sequence)
            event = prepHelper.fileEventMap[sequence]
            if len(prepXmlElems) == 0:
                processedSequences.append(sequence)
            else:
                sequenceJoin = sequence
                while sequenceJoin[0] == '/':
                    sequenceJoin = sequenceJoin[1:]
                processedSequence = os.path.join(tempDir, sequenceJoin)
                processedSequences.append(processedSequence)
                logger.info("Adding child batch_preprocessor target")
                assert sequence != processedSequence
                self.addChildTarget(BatchPreprocessor(self.cactusWorkflowArguments, event, prepXmlElems, 
                                                      sequence, processedSequence, 
                                                      self.getOptionalPhaseAttrib("memory", typeFn=int, default=sys.maxint),
                                                      self.getOptionalPhaseAttrib("cpu", typeFn=int, default=sys.maxint),
                                                      0))
            self.cactusWorkflowArguments.sequences = processedSequences

        # we circumvent makeFollowOnPhaseTarget() interface for this job.
        setupTarget = CactusSetupPhase(cactusWorkflowArguments=self.cactusWorkflowArguments,
                                       phaseName='setup', topFlowerName=self.topFlowerName,
                                       index=0)

        exp = ExperimentWrapper(self.cactusWorkflowArguments.experimentNode)
        if exp.getDbType() == "kyoto_tycoon":
            logger.info("Created ktserver pattern target cactus_setup")
            addKtserverDependentChild(self, setupTarget, isSecondary = False)
        else:
            logger.info("Created follow-on target cactus_setup")
            self.setFollowOnTarget(setupTarget)

############################################################
############################################################
############################################################
##The setup phase.
############################################################
############################################################
############################################################
        
class CactusSetupPhase(CactusPhasesTarget):   
    def run(self):
        runCactusSetup(cactusDiskDatabaseString=self.cactusWorkflowArguments.cactusDiskDatabaseString, 
                       sequences=self.cactusWorkflowArguments.sequences, 
                       newickTreeString=self.cactusWorkflowArguments.speciesTree, 
                       outgroupEvents=self.cactusWorkflowArguments.outgroupEventNames)
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

def getLongestPath(node, distance=0.0):
    """Identify the longest path from the mrca of the leaves of the species tree.
    """
    i, j = distance, distance
    if node.left != None:
        i = getLongestPath(node.left, node.left.distance) + distance
    if node.right != None:  
        j = getLongestPath(node.right, node.right.distance) + distance
    return max(i, j)

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
        if self.getOptionalPhaseAttrib("filterByIdentity", bool, False): #Do the identity filtering
            longestPath = getLongestPath(newickTreeParser(self.cactusWorkflowArguments.speciesTree))
            adjustedPath = float(self.phaseNode.attrib["identityRatio"]) * longestPath + \
            float(self.phaseNode.attrib["minimumDistance"])
            identity = str(100 - int(100 * inverseJukesCantor(adjustedPath)))
            logger.info("The blast stage will filter by identity, the calculated minimum identity is %s from a longest path of %s and an adjusted path of %s" % (identity, longestPath, adjustedPath))
            assert "IDENTITY" in self.phaseNode.attrib["lastzArguments"]
            self.phaseNode.attrib["lastzArguments"] = self.phaseNode.attrib["lastzArguments"].replace("IDENTITY", identity)
        if self.getPhaseIndex() == 0 and self.cactusWorkflowArguments.constraintsFile != None: #Setup the constraints arg
            newConstraintsFile = os.path.join(self.getGlobalTempDir(), "constraints.cig")
            runCactusConvertAlignmentToCactus(self.cactusWorkflowArguments.cactusDiskDatabaseString,
                                              self.cactusWorkflowArguments.constraintsFile, newConstraintsFile)
            self.phaseNode.attrib["constraints"] = newConstraintsFile
            
        if self.getPhaseIndex()+1 < self.getPhaseNumber(): #Check if there is a repeat phase
            self.runPhase(CactusCafRecursion, CactusCafPhase, "caf", index=self.getPhaseIndex()+1)
        else:
            self.runPhase(CactusCafRecursion, CactusBarPhase, "bar")

class CactusCafRecursion(CactusRecursionTarget):
    """This target does the get flowers down pass for the CAF alignment phase.
    """    
    def run(self):
        self.makeRecursiveTargets()
        self.makeExtendingTargets(target=CactusCafWrapper, overlargeTarget=CactusCafWrapperLarge)
        
class CactusCafWrapper(CactusRecursionTarget):
    """Runs cactus_core upon a set of flowers and no alignment file.
    """
    def runCactusCafInWorkflow(self, alignmentFile):
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
                          requiredIngroupFraction=self.getOptionalPhaseAttrib("requiredIngroupFraction", float),
                          requiredOutgroupFraction=self.getOptionalPhaseAttrib("requiredOutgroupFraction", float),
                          requiredAllFraction=self.getOptionalPhaseAttrib("requiredAllFraction", float),
                          singleCopyIngroup=self.getOptionalPhaseAttrib("singleCopyIngroup", bool),
                          singleCopyOutgroup=self.getOptionalPhaseAttrib("singleCopyOutgroup", bool),
                          lastzArguments=self.getOptionalPhaseAttrib("lastzArguments"),
                          minimumSequenceLengthForBlast=self.getOptionalPhaseAttrib("minimumSequenceLengthForBlast", int, 1),
                          maxAdjacencyComponentSizeRatio=self.getOptionalPhaseAttrib("maxAdjacencyComponentSizeRatio", float))
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
        alignmentFile = os.path.join(self.getGlobaTempDir(), "alignments.cigar")
        flowerName = decodeFirstFlowerName(self.flowerNames)
        self.addChildTarget(BlastFlower(self.cactusDiskDatabaseString, 
                                          flowerName, alignmentFile, 
                                          blastOptions=\
                                          BlastOptions(chunkSize=self.getOptionalPhaseAttrib("chunkSize", int),
                                                        overlapSize=self.getOptionalPhaseAttrib("overlapSize", int),
                                                        lastzArguments=self.getOptionalPhaseAttrib("lastzArguments"),
                                                        compressFiles=self.getOptionalPhaseAttrib("compressFiles", bool),
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
        self.makeExtendingTargets(CactusBarWrapper)

class CactusBarWrapper(CactusRecursionTarget):
    """Runs the BAR algorithm implementation.
    """
    def run(self):
        runCactusBar(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                     flowerNames=self.flowerNames, 
                     maximumLength=self.getOptionalPhaseAttrib("bandingLimit", float),
                     spanningTrees=self.getOptionalPhaseAttrib("spanningTrees", int), 
                     gapGamma=self.getOptionalPhaseAttrib( "gapGamma", float), 
                     splitMatrixBiggerThanThis=self.getOptionalPhaseAttrib("splitMatrixBiggerThanThis", int), 
                     anchorMatrixBiggerThanThis=self.getOptionalPhaseAttrib("anchorMatrixBiggerThanThis", int), 
                     repeatMaskMatrixBiggerThanThis=self.getOptionalPhaseAttrib("repeatMaskMatrixBiggerThanThis", int), 
                     diagonalExpansion=self.getOptionalPhaseAttrib("diagonalExpansion"),
                     constraintDiagonalTrim=self.getOptionalPhaseAttrib("constraintDiagonalTrim", int), 
                     minimumBlockDegree=self.getOptionalPhaseAttrib("minimumBlockDegree", int),
                     alignAmbiguityCharacters=self.getOptionalPhaseAttrib("alignAmbiguityCharacters", bool),
                     pruneOutStubAlignments=self.getOptionalPhaseAttrib("pruneOutStrubAlignments", bool),
                     requiredIngroupFraction=self.getOptionalPhaseAttrib("requiredIngroupFraction", float),
                     requiredOutgroupFraction=self.getOptionalPhaseAttrib("requiredOutgroupFraction", float),
                     requiredAllFraction=self.getOptionalPhaseAttrib("requiredAllFraction", float))
        
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
                            maxNumberOfChains=self.getOptionalPhaseAttrib("maxNumberOfChains", int))

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
        self.phaseNode.attrib["secondaryDatabaseString"] = self.cactusWorkflowArguments.secondaryDatabaseString
        self.runPhase(CactusReferenceRecursion, CactusSetReferenceCoordinatesDownPhase, "reference", 
                      doRecursion=self.getOptionalPhaseAttrib("buildReference", bool, False),
                      launchSecondaryKtForRecursiveTarget = True)
        
class CactusReferenceRecursion(CactusRecursionTarget):
    """This target creates the wrappers to run the reference problem algorithm, the follow on target then recurses down.
    """
    def run(self):
        self.makeWrapperTargets(CactusReferenceWrapper)
        self.makeFollowOnRecursiveTarget(CactusReferenceRecursion2)
        
class CactusReferenceWrapper(CactusRecursionTarget):
    """Actually run the reference code.
    """
    def run(self):
        runCactusReference(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                       flowerNames=self.flowerNames, 
                       matchingAlgorithm=self.getOptionalPhaseAttrib("matchingAlgorithm"), 
                       permutations=self.getOptionalPhaseAttrib("permutations", int),
                       referenceEventString=self.getOptionalPhaseAttrib("reference"), 
                       useSimulatedAnnealing=self.getOptionalPhaseAttrib("useSimulatedAnnealing", bool),
                       theta=self.getOptionalPhaseAttrib("theta", float),
                       maxWalkForCalculatingZ=self.getOptionalPhaseAttrib("maxWalkForCalculatingZ", int))
        
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
    def getMafFile(self):
        return self.cactusWorkflowArguments.experimentNode.find("hal").attrib["mafPath"]
        
    def getHalFile(self):
        return self.cactusWorkflowArguments.experimentNode.find("hal").attrib["halPath"]
    
    def raiseHal(self, makeMaf):
        self.setupSecondaryDatabase()
        self.phaseNode.attrib["secondaryDatabaseString"] = self.cactusWorkflowArguments.secondaryDatabaseString
        if makeMaf:
            self.phaseNode.attrib["outputFile"]=self.getMafFile()
            self.phaseNode.attrib["makeMaf"]="1"
            self.makeFollowOnPhaseTarget(CactusHalGeneratorPhaseCleanup, "hal")
        else:
            self.phaseNode.attrib["outputFile"]=self.getHalFile()
            self.phaseNode.attrib["makeMaf"]="0"
            self.makeFollowOnPhaseTarget(CactusHalGeneratorPhase2, "hal")
        self.makeRecursiveChildTarget(CactusHalGeneratorRecursion, launchSecondaryKtForRecursiveTarget=True)
    
    def run(self):
        referenceNode = findRequiredNode(self.cactusWorkflowArguments.configNode, "reference")
        if referenceNode.attrib.has_key("reference"):
            self.phaseNode.attrib["reference"] = referenceNode.attrib["reference"]
        if self.getOptionalPhaseAttrib("buildFasta", bool, default=False):
            self.phaseNode.attrib["fastaPath"] = self.cactusWorkflowArguments.experimentNode.find("hal").attrib["fastaPath"]
            self.makeRecursiveChildTarget(CactusFastaGenerator)
        if self.getOptionalPhaseAttrib("buildHal", bool, default=False):
            self.raiseHal(makeMaf=False)
        elif self.getOptionalPhaseAttrib("buildMaf", bool, default=False):
            self.raiseHal(makeMaf=True)

class CactusFastaGenerator(CactusRecursionTarget):
    def run(self):
        runCactusFastaGenerator(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                                    flowerName=decodeFirstFlowerName(self.flowerNames),
                                    outputFile=self.getOptionalPhaseAttrib("fastaPath"),
                                    referenceEventString=self.getOptionalPhaseAttrib("reference"))
            
class CactusHalGeneratorPhase2(CactusHalGeneratorPhase):
    def run(self): #This allows us to generate a maf
        if self.getOptionalPhaseAttrib("buildHal", bool, default=False):
            self.cleanupSecondaryDatabase()
        if self.getOptionalPhaseAttrib("buildMaf", bool, default=False):
            self.raiseHal(makeMaf=True)

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
                              self.getOptionalPhaseAttrib("showOnlySubstitutionsWithRespectToReference", bool),
                              makeMaf=self.getOptionalPhaseAttrib("makeMaf", bool))

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
        self.experimentFile = options.experimentFile
        self.experimentNode = ET.parse(self.experimentFile).getroot()
        #Get the database string
        self.cactusDiskDatabaseString = ET.tostring(self.experimentNode.find("cactus_disk").find("st_kv_database_conf"))
        #Get the species tree
        self.speciesTree = self.experimentNode.attrib["species_tree"]
        #Get the sequences
        self.sequences = self.experimentNode.attrib["sequences"].split()
        #Get any list of 'required species' for the blocks of the cactus.
        self.outgroupEventNames = getOptionalAttrib(self.experimentNode, "outgroup_events")
        #Constraints
        self.constraintsFile = getOptionalAttrib(self.experimentNode, "constraints")
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
            
        #The config options
        configFile = self.experimentNode.attrib["config"]
        if configFile == "default":
            configFile = os.path.join(cactusRootPath(), "pipeline", "cactus_workflow_config.xml")
        elif configFile == "defaultProgressive":
            configFile = os.path.join(cactusRootPath(), "progressive", "cactus_progressive_workflow_config.xml")
        else:
            logger.info("Using user specified config file: %s", configFile)
        self.configNode = ET.parse(configFile).getroot()
        if options.buildAvgs:
            findRequiredNode(self.configNode, "avg").attrib["buildAvgs"] = "1"
        if options.buildReference:
            findRequiredNode(self.configNode, "reference").attrib["buildReference"] = "1"
        if options.buildHal:
            findRequiredNode(self.configNode, "hal").attrib["buildHal"] = "1"
        if options.buildMaf:
            findRequiredNode(self.configNode, "hal").attrib["buildMaf"] = "1"
        if options.buildFasta:
            findRequiredNode(self.configNode, "hal").attrib["buildFasta"] = "1"

def addCactusWorkflowOptions(parser):
    parser.add_option("--experiment", dest="experimentFile", 
                      help="The file containing a link to the experiment parameters")
    
    parser.add_option("--skipAlignments", dest="skipAlignments", action="store_true",
                      help="Skip building alignments", default=False)
    
    parser.add_option("--buildAvgs", dest="buildAvgs", action="store_true",
                      help="Build trees", default=False)
    
    parser.add_option("--buildReference", dest="buildReference", action="store_true",
                      help="Creates a reference ordering for the flowers", default=False)
    
    parser.add_option("--buildHal", dest="buildHal", action="store_true",
                      help="Build a hal file", default=False)
    
    parser.add_option("--buildMaf", dest="buildMaf", action="store_true",
                      help="Build a maf file", default=False)
    
    parser.add_option("--buildFasta", dest="buildFasta", action="store_true",
                      help="Build a fasta file of the input sequences (and reference sequence, used with hal output)", default=False)
    
def main():
    ##########################################
    #Construct the arguments.
    ##########################################
    
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)
    addCactusWorkflowOptions(parser)
        
    options, args = parser.parse_args()
    setLoggingFromOptions(options)
    
    if len(args) != 0:
        raise RuntimeError("Unrecognised input arguments: %s" % " ".join(args))

    cactusWorkflowArguments = CactusWorkflowArguments(options)
    if options.skipAlignments: #Don't dp caf and bar, jump right to an existing cactus structure
        #cactusWorkflowArguments, phaseName, topFlowerName, index=0):
        firstTarget = CactusNormalPhase(cactusWorkflowArguments=cactusWorkflowArguments, phaseName="normal")
    else:
        firstTarget = CactusSetupPhase(cactusWorkflowArguments=cactusWorkflowArguments, phaseName="setup")
    Stack(firstTarget).startJobTree(options)

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.pipeline.cactus_workflow import *
    _test()
    main()
