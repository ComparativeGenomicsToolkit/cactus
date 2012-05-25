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
from optparse import OptionParser

from sonLib.bioio import getTempFile
from sonLib.bioio import newickTreeParser

from sonLib.bioio import logger
from sonLib.bioio import setLoggingFromOptions
from sonLib.bioio import getTempDirectory

from cactus.shared.common import cactusRootPath
  
from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack 

from cactus.shared.common import runCactusSetup
from cactus.shared.common import runCactusCaf
from cactus.shared.common import runCactusGetFlowers
from cactus.shared.common import runCactusExtendFlowers
from cactus.shared.common import encodeFlowerNames
from cactus.shared.common import decodeFlowerNames
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

from cactus.blastAlignment.cactus_aligner import MakeSequences
from cactus.blastAlignment.cactus_batch import MakeBlastOptions
from cactus.blastAlignment.cactus_batch import makeBlastFromOptions

from cactus.preprocessor.cactus_preprocessor import BatchPreprocessor
from cactus.preprocessor.cactus_preprocessor import PreprocessorHelper

############################################################
############################################################
############################################################
##The preprocessor phase.
############################################################
############################################################
############################################################

class CactusPreprocessorPhase(Target):
    def __init__(self, options, sequences):
        Target.__init__(self, time=0.0002)
        self.options = options 
        self.sequences = sequences

    def run(self):
        self.logToMaster("Starting preprocessor phase target at %s seconds" % time.time())
        tempDir = getTempDirectory(self.getGlobalTempDir())
        prepHelper = PreprocessorHelper(self.options, self.sequences)
        processedSequences = []
        for sequence in self.sequences:
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
                self.addChildTarget(BatchPreprocessor(self.options, event, prepXmlElems, sequence, processedSequence, 0))
        self.setFollowOnTarget(CactusSetupPhase(self.options, processedSequences))
        logger.info("Created followOn target cactus_setup job, and follow on down pass job")

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
    nodes = configNode.findAll(nodeName)
    if nodes == None:
        raise RuntimeError("Could not find any nodes with name %s in %s node" % (nodeName, configNode))
    if index >= len(nodes):
        raise RuntimeError("Could not find a node with name %s and index %i in %s node" % (nodeName, index, configNode))
    return nodes[index]

def extractNode(node):
    """Make an XML node free of its parent subtree
    """
    return ET.fromstring(ET.tostring(node))

class CactusTarget(Target):
    """Base target for all cactus workflow targets.
    """
    def __init__(self, phaseNode, overlarge=False):
        self.phaseNode = phaseNode
        self.overlarge = overlarge
        className = str(self.__class__).split(".")[-1]
        assert className != ''
        self.targetNode = self.phaseNode.find(className)
        if overlarge:
            Target.__init__(self, memory=self.getOptionalTargetAttrib("overlargeMemory", sys.maxint), 
                            cpu=self.getOptionalTargetAttrib("overlargeCpu", sys.maxint))
        else:
            Target.__init__(self, memory=self.getOptionalTargetAttrib("memory", sys.maxint), 
                            cpu=self.getOptionalTargetAttrib("cpu", sys.maxint))
    
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
    def __init__(self, options, phaseName, topFlowerName, index=0):
        phaseNode = findRequiredNode(options.configNode, phaseName, index)
        CactusTarget.__init__(self, phaseNode=phaseNode)
        self.index = index
        self.options = options
        self.topFlowerName = topFlowerName
    
    def makeRecursiveChildTarget(self, target):
        self.addChildTarget(target(phaseNode=extractNode(self.phaseNode), 
                                   cactusDiskDatabaseString=self.options.cactusDiskDatabaseString, 
                                   flowerNames=encodeFlowerNames((self.topFlowerName,)), overlarge=True))
    
    def makeFollowOnPhaseTarget(self, target, phaseName, index=0):
        self.setFollowOnTarget(target(options=self.options, phaseName=phaseName, topFlowerName=self.topFlowerName, index=index))
        
    def runPhase(self, recursiveTarget, nextPhaseTarget, nextPhaseName, doRecursion=True):
        self.logToMaster("Starting %s phase target at %s seconds" % (self.phaseNode.name, time.time()))
        if doRecursion:
            self.makeRecursiveChildTarget(recursiveTarget)
        self.makeFollowOnPhaseTarget(nextPhaseTarget, nextPhaseName)
        
    def getPhaseIndex(self):
        return self.index
    
    def getPhaseNumber(self):
        return len(self.options.configNode.findAll(self.phaseNode.tag))

class CactusRecursionTarget(CactusTarget):
    """Base recursive target for traversals up and down the cactus tree.
    """
    MAX_SEQUENCE_SIZE_OF_FLOWER_GROUPING = 1000000
    def __init__(self, phaseNode, cactusDiskDatabaseString, flowerNames, overlarge=False):
        CactusTarget.__init__(self, phaseNode=phaseNode, overlarge=overlarge)
        self.cactusDiskDatabaseString = cactusDiskDatabaseString
        self.flowerNames = flowerNames  
        
    def makeFollowOnRecursiveTarget(target):
        """Sets the followon to the given recursive target
        """
        self.setFollowOnTarget(target(phaseNode=phaseNode, 
                                   cactusDiskDatabaseString=self.options.cactusDiskDatabaseString, 
                                   flowerNames=encodeFlowerNames((self.topFlowerName,)), overlarge=self.overlarge))
        
    def makeChildTargets(flowersAndSizes, target, overlargeTarget=None, phaseNode=None):
        """Make a set of child targets for a given set of flowers and chosen child target
        """
        if overlargeTarget == None:
            overlargeTarget = target
        if phaseNode == None:
            phaseNode = self.phaseNode
        for overlarge, flowerNames in flowersAndSizes:
            if overlarge: #Make sure large flowers are on there own, in their own job
                flowerStatsString = runCactusFlowerStats(cactusDiskDatabaseString=self.cactusDiskDatabaseString, decodeFlowerNames(flowerNames)[0])
                self.logToMaster("Adding an oversize flower for target class %s and stats %s" \
                                         % (overlargeTarget.__class__, flowerStatsString))
                self.addChildTarget(overlargeTarget(cactusDiskDatabaseString=self.cactusDiskDatabaseString, phaseNode=phaseNode, 
                                                    flowerNames=flowerNames, overlarge=1)) #This ensures overlarge flowers, 
            else:
                self.addChildTarget(target(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                                           phaseNode=phaseNode, flowerNames=flowerNames))
        
    def makeRecursiveTargets(self, phaseNode=None):
        """Make a set of child targets for a given set of parent flowers.
        """
        self.makeChildTargets(flowersAndSizes=runCactusGetFlowers(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                                                                  flowerNames=self.flowerNames, 
                                                                  minSequenceSizeOfFlower=self.getOptionalTargetAttrib("minFlowerSize", int, 0),
                                                                  maxSequenceSizeOfFlowerGrouping=self.getOptionalTargetAttrib("maxFlowerGroupSize", 
                                                                                                                int, default=CactusRecursionTarget.MAX_SEQUENCE_SIZE_OF_FLOWER_GROUPING)), 
                              target=self.__class__, phaseNode=phaseNode) 
    
    def makeExtendingTargets(self, target, overlargeTarget=None, phaseNode=None):
        """Make set of child targets that extend the current cactus tree.
        """
        self.makeChildTargets(flowersAndSizes=runCactusExtendFlowers(self.cactusDiskDatabaseString, self.flowerNames, 
                                           minSequenceSizeOfFlower=self.getOptionalTargetAttrib("minFlowerSize", int, 1), 
                                           maxSequenceSizeOfFlowerGrouping=self.getOptionalTargetAttrib("maxFlowerGroupSize", int, 
                                           default=CactusRecursionTarget.MAX_SEQUENCE_SIZE_OF_FLOWER_GROUPING)), 
                              parentTarget=self, target=target,
                              overlargeTarget=overlargeTarget, phaseNode=phaseNode)
    
    def makeWrapperTargets(self, target, overlargeTarget=None, phaseNode=None):
        """Takes the list of flowers for a recursive target and splits them up to fit the given wrapper target(s).
        """
        self.makeChildTargets(flowersAndSizes=runCactusSplitFlowersList(self.cactusDiskDatabaseString, self.flowerNames, 
                                           maxSequenceSizeOfFlowerGrouping=self.getOptionalTargetAttrib("maxFlowerGroupSize", int, 
                                           default=CactusRecursionTarget.MAX_SEQUENCE_SIZE_OF_FLOWER_GROUPING)), 
                              parentTarget=self, target=target,
                              overlargeTarget=overlargeTarget, phaseNode=phaseNode)

############################################################
############################################################
############################################################
##The setup phase.
############################################################
############################################################
############################################################
        
class CactusSetupPhase(CactusPhasesTarget):   
    def run(self):
        self.runPhase(CactusSetupWrapper, CactusCafPhase, "caf")
        
class CactusSetupWrapper(CactusRecursionTarget):
    def run(self):
        runCactusSetup(cactusDiskDatabaseString=self.cactusDiskDatabaseString, sequences=self.options.sequences, 
                       speciesTree=self.options.speciesTree, outgroupEvents=self.options.outgroupEventNames)

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
    """Identify the longest path from root to leaves of the species tree
    and add the min-distance.
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
    
    >>> inverseJukesCantor(0.5)
    0.36493716072555599
    >>> inverseJukesCantor(1.0)
    0.55230214641320496
    >>> inverseJukesCantor(10.0)
    0.74999878530240571
    >>> inverseJukesCantor(100000.0)
    0.75
    """
    assert d >= 0.0
    return 0.75 * (1 - math.exp(-d * 4.0/3.0))
    
class CactusCafPhase(CactusPhasesTarget):      
    def run(self):
        if self.getOptionalPhaseAttrib("filterByIdentity", bool, False): #Do the identity filtering
            longestPath = getLongestPath(newickTreeParser(self.options.speciesTree))
            adjustedPath = float(self.phaseNode.attrib["identityRatio"]) * longestPath + \
            float(self.phaseNode.attrib["minimumDistance"])
            identity = str(100 - int(100 * inverseJukesCantor(adjustedPath)))
            logger.info("The blast stage will filter by identity, the calculated minimum identity is %s from a longest path of %s and an adjusted path of %s" % (identity, longestPath, adjustedPath))
            assert "IDENTITY" in self.phaseNode.attrib["lastzArguments"]
            self.phaseNode.attrib["lastzArguments"] = self.phaseNode.attrib["lastzArguments"].replace("IDENTITY", identity)
        if self.phaseIndex() == 0 and "constraints" in self.options.experimentFile.attrib: #Setup the constraints arg
            newConstraintsFile = os.path.join(self.getGlobalTempDir(), "constraints.cig")
            runCactusConvertAlignmentToCactus(self.options.cactusDiskDatabaseString,
                                              self.options.experimentFile.attrib["constraints"], newConstraintsFile)
            self.phaseNode.attrib["constraints"] = newConstraintsFile
        if self.getPhaseIndex() < self.getPhaseNumber("caf"): #Check if there is a repeat phase
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
        self.runCactusCafInWorkflow(None)
       
class CactusCafWrapperLarge(CactusRecursionTarget):
    """Runs blast on the given flower and passes the resulting alignment to cactus core.
    """
    def run(self):
        logger.info("Starting the cactus aligner target")
        #Generate a temporary file to hold the alignments
        alignmentFile = getTempFile(".fa", self.getGlobalTempDir())
        logger.info("Got an alignments file")
        #Now make the child aligner target
        flowerNames = decodeFlowerNames(self.flowerNames)
        assert len(flowerNames) == 1
        blastOptions = \
        makeBlastFromOptions(MakeBlastOptions(chunkSize=self.getOptionalPhaseAttrib("chunkSize", int),
                                              overlapSize=self.getOptionalPhaseAttrib("overlapSize", int),
                                              lastzArguments=self.getOptionalPhaseAttrib("lastzArguments"),
                                              chunksPerJob=self.getOptionalPhaseAttrib("chunksPerJob", int),
                                              compressFiles=self.getOptionalPhaseAttrib("compressFiles", bool)))
        self.addChildTarget(MakeSequences(self.cactusDiskDatabaseString, 
                                          flowerNames[0], alignmentFile, blastOptions=blastOptions,
                                          minimumSequenceLength=self.getOptionalPhaseAttrib("minimumSequenceLengthForBlast", int, 1)))
        logger.info("Created the cactus_aligner child target")
        #Now setup a call to cactus core wrapper as a follow on
        self.phaseNode.attrib["alignments"] = alignmentFile
        self.makeFollowOnRecursiveTarget(CactusCafWrapper2)
        logger.info("Setup the follow on cactus_core target")
        
class CactusCafWrapperLarge2(CactusCafWrapper):
    """Runs cactus_core upon a one flower and one alignment file.
    """
    def run(self):
        self.runCactusCafInWorkflow(self.phaseNode.attrib["alignments"])
        
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
            self.phaseNode["iterations"] = str(normalisationIterations-1)
            self.runPhase(CactusNormalRecursion, CactusNormalPhase, "normal")
        else:
            self.makeFollowOnPhaseTarget(CactusAVGPhase, "avg")
     
class CactusNormalRecursion(CactusRecursionTarget):
    """This target does the down pass for the normal phase.
    """
    def run(self):
        self.makeRecursiveTargets()
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
        self.makeRecursiveTargets()
        self.makeWrapperTargets(CactusAVGWrapper)

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
        self.runPhase(CactusReferenceRecursion, CactusSetReferenceCoordinatesDownPhase, "reference", 
                      doRecursive=self.getOptionalPhaseAttrib("buildReference", bool, False))
        
class CactusReferenceRecursion(CactusRecursionTarget):
    """This target creates the wrappers to run the reference problem algorithm, the follow on target then recurses down.
    """
    def run(self):
        self.makeWrapperTargets(CactusReferenceWrapper)
        self.makeFollowOnRecursiveTarget(CactusReferenceRecursion2)
        
class CactusReferenceWrapper(CactusRecursionTarget):
    """Actually run the reference code.
    """
    runCactusReference(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                       flowerNames=self.flowerNames, 
                       matchingAlgorithm=self.getOptionalPhaseAttrib("matchingAlgorithm"), 
                       permutations=self.getOptionalPhaseAttrib("permutations", int),
                       referenceEventString=self.getOptionalPhaseAttrib("reference"), 
                       useSimulatedAnnealing=self.getOptionalPhaseAttrib("useSimulatedAnnealing", bool),
                       theta=self.getOptionalPhaseAttrib("theta", float),
                       maxNumberOfChainsBeforeSwitchingToFast=self.getOptionalPhaseAttrib("maxNumberOfChainsBeforeSwitchingToFast", int))
        
class CactusReferenceRecursion2(CactusRecursionTarget):
    def run(self):
        self.makeRecursiveTargets()
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
                                         flowerNames=self.flowerNames,
                                         referenceEventString=self.getOptionalAttrib("reference"), 
                                         outgroupEventString=self.getOptionalAttrib("outgroup"), 
                                         bottomUpPhase=True)
        
class CactusSetReferenceCoordinatesDownPhase(CactusPhasesTarget):
    """This is the second part of the reference coordinate up pass.
    """
    def run(self):
        self.runPhase(CactusSetReferenceCoordinatesDownRecursion, CactusCheckPhase, "check", doRecursive=self.getOptionalAttrib("buildReference", bool, False))
        
class CactusSetReferenceCoordinatesDownRecursion(CactusRecursionTarget):
    """Does the down pass for filling Fills in the coordinates, once a reference is added.
    """        
    def run(self):
        self.makeRecursiveTargets()
        self.makeFollowOnRecursiveTarget(CactusSetReferenceCoordinatesRecursion2)

class CactusSetReferenceCoordinatesRecursion2(CactusRecursionTarget):
    def run(self):
        self.makeWrapperTargets(CactusSetReferenceCoordinatesDownWrapper)

class CactusSetReferenceCoordinatesDownWrapper(CactusRecursionTarget):
    """Does the down pass for filling Fills in the coordinates, once a reference is added.
    """        
    def run(self):
        runCactusAddReferenceCoordinates(cactusDiskDatabaseString=self.cactusDiskDatabaseString, flowerNames=self.flowerNames,
                                         referenceEventString=self.getOptionalPhaseAttrib("reference"),
                                         outgroupEventString=self.getOptionalPhaseAttrib("outgroup"), 
                                         bottomUpPhase=False)

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
        self.runPhase(CactusCheckRecursion, CactusHalGeneratorPhase, "hal", doRecursive=self.getOptionalAttrib("runCheck", bool, False))
        
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
        runCactusCheck(self.cactusDiskDatabaseString, self.flowerNames, self.getOptionalAttrib("checkNormalised", bool, False))

############################################################
############################################################
############################################################
#Hal generation
############################################################
############################################################
############################################################

class CactusHalGeneratorPhase(CactusPhasesTarget):
    def run(self):
        self.logToMaster("Starting the hal generation phase at %s seconds" % time.time())
        if self.getOptionalPhaseAttrib("buildHal", bool, default=False):
            referenceNode = getRequiredNode(self.options.configNode, "reference")
            if referenceNode.attrib.has_key("reference"):
                self.phaseNode.attrib["reference"] = referenceNode.attrib["reference"]
            self.makeRecursiveChildTarget(CactusHalGeneratorRecursion)
            self.phaseNode.attrib["outputFile"] = self.options.experimentFile.find("hal").attrib["path"]
            self.makeRecursiveChildTarget(CactusHalGeneratorRecursion)

class CactusHalGeneratorRecursion(CactusRecursionTarget):
    """Generate the hal file by merging indexed hal files from the children.
    """ 
    def run(self):
        i = extractNode(self.phaseNode)
        i.attrib["parentDir"] = self.getGlobalTempDir()
        if "outputFile" in i.attrib:
            i.attrib.pop("outputFile")
        self.makeRecursiveTargets(phaseNode=i)
        self.makeFollowOnRecursiveTarget(CactusHalGeneratorUpRunnable)
        
class CactusHalGeneratorUpRunnable(CactusHalGeneratorUp):
    """Does the up pass for filling in the coordinates, once a reference is added.
    """ 
    def run(self):
        runCactusHalGenerator(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                              flowerNames=self.flowerNames,
                              referenceEventString=self.getOptionalPhaseAttrib("reference"), #self.configNode.attrib["reference"], #self.getOptionalPhaseAttrib("reference"), 
                              childDir=self.getGlobalTempDir(), 
                              parentDir=self.getOptionalPhaseAttrib("parentDir"),
                              outputFile=self.getOptionalPhaseAttrib("outputFile"),
                              showOnlySubstitutionsWithRespectToReference=\
                              self.getOptionalPhaseAttrib("showOnlySubstitutionsWithRespectToReference", bool),
                              makeMaf=self.getOptionalPhaseAttrib("makeMaf", bool))

# add stuff to the options object 
# (code extracted from the main() method so it can be reused by progressive)
def expandWorkflowOptions(options, experimentFile = None):
    if experimentFile is not None:
        options.experimentFile = experimentFile
    else:
        options.experimentFile = ET.parse(options.experimentFile).getroot()
    #Get the database string
    options.cactusDiskDatabaseString = ET.tostring(options.experimentFile.find("cactus_disk").find("st_kv_database_conf"))
    #Get the species tree
    options.speciesTree = options.experimentFile.attrib["species_tree"]
    #Parse the config file which contains all the program options
    if options.experimentFile.attrib["config"] == "default":
        options.experimentFile.attrib["config"] = os.path.join(cactusRootPath(), "pipeline", "cactus_workflow_config.xml")
    else:
        logger.info("Using user specified config file")
    #Get the config file for the experiment
    options.config = ET.parse(options.experimentFile.attrib["config"]).getroot()
    #Get any list of 'required species' for the blocks of the cactus.
    options.outgroupEventNames = getOptionalAttrib(options.experimentFile, "outgroup_events")
    logger.info("Parsed the XML options file")
    
def main():
    ##########################################
    #Construct the arguments.
    ##########################################
    
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)
    
    parser.add_option("--experiment", dest="experimentFile", help="The file containing a link to the experiment parameters")
    
    parser.add_option("--setupAndBuildAlignments", dest="setupAndBuildAlignments", action="store_true",
                      help="Setup and build alignments then normalise the resulting structure", default=False)
    
    parser.add_option("--buildTrees", dest="buildTrees", action="store_true",
                      help="Build trees", default=False) 
    
    parser.add_option("--buildFaces", dest="buildFaces", action="store_true",
                      help="Build adjacencies", default=False)
    
    parser.add_option("--buildReference", dest="buildReference", action="store_true",
                      help="Creates a reference ordering for the flowers", default=False)
    
    parser.add_option("--buildHal", dest="buildHal", action="store_true",
                      help="Build a hal file", default=False)
    
    options, args = parser.parse_args()
    setLoggingFromOptions(options)
    
    if len(args) != 0:
        raise RuntimeError("Unrecognised input arguments: %s" % " ".join(args))

    # process the options
    expandWorkflowOptions(options)
    #Get the sequences
    sequences = options.experimentFile.attrib["sequences"].split()

    if options.setupAndBuildAlignments:
        baseTarget = CactusSetupPhase(options, sequences)
        logger.info("Going to create alignments and define the cactus tree")
    elif options.buildTrees:
        baseTarget = CactusPhylogenyPhase('0', options)
        logger.info("Starting from phylogeny phase")
    elif options.buildReference:
        baseTarget = CactusReferencePhase('0', options)
        logger.info("Starting from reference phase")
    elif options.buildFaces:
        baseTarget = CactusFacesPhase('0', options)
        logger.info("Starting from faces phase")
    else:
        logger.info("Nothing to do!")
        return
    
    Stack(baseTarget).startJobTree(options)
    logger.info("Done with job tree")

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.pipeline.cactus_workflow import *
    _test()
    main()
