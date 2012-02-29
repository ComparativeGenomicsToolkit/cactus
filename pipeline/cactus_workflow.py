#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Script strings together all the components to make the basic pipeline for reconstruction.

The script uses the the jobTree.scriptTree target framework so structure all the related wrappers.

There are four high level wrappers, a SetupPhase, DownPassPhase, UpPassPhase, VerificationPahse. 

In the setup phase the system sets up the files needed for the reconstruction problem.

In the down pass phase alignments and trees are built.

In the up pass phase the adjacencies are added.

In the verification phase the reconstruction tree is checked against the expected spec.

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
from cactus.shared.common import runCactusCore
from cactus.shared.common import runCactusGetFlowers
from cactus.shared.common import runCactusExtendFlowers
from cactus.shared.common import runCactusPhylogeny
from cactus.shared.common import runCactusAdjacencies
from cactus.shared.common import runCactusBaseAligner
from cactus.shared.common import runCactusMakeNormal 
from cactus.shared.common import runCactusReference
from cactus.shared.common import runCactusAddReferenceCoordinates
from cactus.shared.common import runCactusCheck
from cactus.shared.common import runCactusRecursiveMafGenerator
from cactus.shared.common import runCactusFlowerStats

from cactus.blastAlignment.cactus_aligner import MakeSequences
from cactus.blastAlignment.cactus_batch import MakeBlastOptions
from cactus.blastAlignment.cactus_batch import makeBlastFromOptions

from cactus.preprocessor.cactus_preprocessor import BatchPreprocessor
from cactus.preprocessor.cactus_preprocessor import PreprocessorHelper

############################################################
############################################################
############################################################
##The setup phase.
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

def getOptionalAttrib(node, attribName, typeFn=None, default=None):
    """Get an optional attrib, or None, if not set.
    """
    if node.attrib.has_key(attribName):
        if typeFn != None:
            if typeFn == bool:
                return typeFn(int(node.attrib[attribName]))
            return typeFn(node.attrib[attribName])
        return node.attrib[attribName]
    return default
        
class CactusSetupPhase(Target):
    def __init__(self, options, sequences):
        Target.__init__(self, time=0.0002)
        self.options = options
        self.sequences = sequences
            
    def modifyConfig(self):
        #Add the identity clause into the blast strings
        alignmentNode = self.options.config.find("alignment")
        for iterationNode in alignmentNode.find("iterations").findall("iteration"):
            blastNode = iterationNode.find("blast")
            if blastNode != None:
                assert iterationNode.attrib["type"] == "blast"
                if getOptionalAttrib(blastNode, "filterByIdentity", bool, False):
                    longestPath = getLongestPath(newickTreeParser(self.options.speciesTree))
                    adjustedPath = float(blastNode.attrib["identityRatio"]) * longestPath + \
                    float(blastNode.attrib["minimumDistance"])
                    identity = str(100 - int(100 * inverseJukesCantor(adjustedPath)))
                    logger.info("The blast stage will filter by identity, the calculated minimum identity is %s from a longest path of %s and an adjusted path of %s" % (identity, longestPath, adjustedPath))
                    assert "IDENTITY" in blastNode.attrib["lastzArguments"]
                    blastNode.attrib["lastzArguments"] = blastNode.attrib["lastzArguments"].replace("IDENTITY", identity)
                    
    def run(self):
        self.logToMaster("Starting setup phase target at %s seconds" % time.time())
        #Modify the config options
        self.modifyConfig()
        #Make the child setup job.
        self.addChildTarget(CactusPreprocessorPhase(self.options, self.sequences))
        #initialise the down pass as the follow on.. using special '0'
        self.setFollowOnTarget(CactusAlignmentPhase(self.options, 0))
        logger.info("Created child target preprocessor job, and follow on down pass job")

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

        self.setFollowOnTarget(CactusSetupWrapper(self.options, processedSequences))
        logger.info("Created followOn target cactus_setup job, and follow on down pass job")
        
class CactusSetupWrapper(Target):
    def __init__(self, options, sequences):
        Target.__init__(self)
        self.options = options
        self.sequences = sequences
        
    def run(self):
        logger.info("Starting cactus setup target")
        runCactusSetup(self.options.cactusDiskDatabaseString, self.sequences, 
                       self.options.speciesTree, outgroupEvents=self.options.outgroupEventNames)
        logger.info("Finished the setup phase target")

############################################################
############################################################
############################################################
#The alignment phases, split into the Caf and Bar phases.
############################################################
############################################################
############################################################

def extractNode(node):
    """Make an XML node free of its parent subtree
    """
    return ET.fromstring(ET.tostring(node))

class CactusPhasesTarget(Target):
    """Base target for each phase.
    """
    def __init__(self, options, flowerName, iteration=0):
        Target.__init__(self, time=0.0002)
        self.options = options
        self.flowerName = flowerName
        self.iteration = iteration
    
class CactusAlignmentPhase(CactusPhasesTarget):        
    def run(self):
        self.logToMaster("Starting the alignment phase for iteration %i at %i seconds" % (self.iteration, time.time()))
        iterations = self.options.config.find("alignment").find("iterations").findall("iteration")
        if self.iteration < len(iterations):
            configNode = extractNode(iterations[self.iteration])
            if configNode.attrib["type"] == "blast":
                self.addChildTarget(CactusCafDown(self.options.cactusDiskDatabaseString, configNode, 
                                                  [ self.flowerName, 1 ]))
            else:
                assert configNode.attrib["type"] == "base"
                self.addChildTarget(CactusBarDown(self.options.cactusDiskDatabaseString, configNode, 
                                                  [ self.flowerName, 1 ]))
            self.setFollowOnTarget(CactusAlignmentPhase(self.options, self.flowerName, self.iteration+1))
        else:
            self.setFollowOnTarget(CactusNormalPhase(self.options, self.flowerName))

############################################################
############################################################
############################################################
#The CAF phase.
#
#Creates the reconstruction structure with blocks
############################################################
############################################################
############################################################

class CactusRecursionTarget(Target):
    """Base recursive target for traversals up and down the cactus tree.
    """
    def __init__(self, cactusDiskDatabaseString, configNode, flowerNames, memory=sys.maxint, cpu=sys.maxint):
        targetStringId = str(self.__class__)[8:-2]
        Target.__init__(self, memory=memory, cpu=cpu)
        self.cactusDiskDatabaseString = cactusDiskDatabaseString
        self.configNode = configNode
        self.flowerNames = flowerNames
  
MAX_SEQUENCE_SIZE_OF_FLOWER_GROUPING=1000000        

def makeTargets(cactusDiskDatabaseString, configNode, flowersAndSizes, 
                parentTarget, target, overlargeTarget=None, 
                memory=sys.maxint, cpu=sys.maxint, 
                overlargeMemory=sys.maxint, overlargeCpu=sys.maxint):
    """Make a set of targets for a given set of flowers.
    """
    flowerNames = []
    totalSequenceSize = 0.0
    maxSequenceSizeOfFlowerGrouping=getOptionalAttrib(configNode, "maxFlowerGroupSize", int, default=MAX_SEQUENCE_SIZE_OF_FLOWER_GROUPING)
    flowerGrouping = []
    if overlargeTarget == None:
        overlargeTarget = target
    for totalFlowerSize, firstFlowerName, flowerNumber in flowersAndSizes:
        if totalFlowerSize > maxSequenceSizeOfFlowerGrouping: #Make sure large flowers are on there own, in their own job
            assert flowerNumber == 1
            flowerStatsString = runCactusFlowerStats(cactusDiskDatabaseString, firstFlowerName)
            parentTarget.logToMaster("Adding an oversize flower: %s on its own, with %s bases for target class %s and stats %s" \
                                     % (firstFlowerName, totalFlowerSize, target, flowerStatsString))
            parentTarget.addChildTarget(overlargeTarget(cactusDiskDatabaseString=cactusDiskDatabaseString, 
                                                        configNode=configNode, 
                                                        flowerNames=[ firstFlowerName, 1 ], memory=overlargeCpu, cpu=overlargeCpu)) #This ensures overlarge flowers, 
            #an in cactus core, get diverted and run on their own.
        else:
            totalSequenceSize += totalFlowerSize
            flowerNames.append(firstFlowerName)
            flowerNames.append(flowerNumber)
            if totalSequenceSize >= maxSequenceSizeOfFlowerGrouping: 
                flowerGrouping.append(flowerNames)
                flowerNames = []
                totalSequenceSize = 0.0
    if len(flowerNames) > 0:
        assert totalSequenceSize < maxSequenceSizeOfFlowerGrouping
        if len(flowerGrouping) > 0: #Avoid small targets if multiple targets exist by adding remaining targets jobs
            k = 0
            l = len(flowerGrouping)
            while len(flowerNames) > 0:
                i = flowerNames.pop()
                j = flowerNames.pop()
                flowerGrouping[k % l].append(j)
                flowerGrouping[k % l].append(i)
                k += 1
        else:
            flowerGrouping.append(flowerNames)
    for flowerNames in flowerGrouping:
        parentTarget.addChildTarget(target(cactusDiskDatabaseString=cactusDiskDatabaseString, 
                                           configNode=configNode, 
                                           flowerNames=flowerNames, 
                                           memory=memory, cpu=cpu))
     
def makeChildTargets(cactusDiskDatabaseString, configNode, flowerNames, target, childTarget):
    """Make a set of child targets for a given set of parent flowers.
    """
    minSequenceSizeOfFlower = getOptionalAttrib(configNode, "minFlowerSize", int, 0) 
    maxSequenceSizeOfFlowerGrouping=getOptionalAttrib(configNode, "maxFlowerGroupSize", int, default=MAX_SEQUENCE_SIZE_OF_FLOWER_GROUPING)
    childFlowers = runCactusGetFlowers(cactusDiskDatabaseString, flowerNames, 
                                       minSequenceSizeOfFlower=minSequenceSizeOfFlower,
                                       maxSequenceSizeOfFlowerGrouping=maxSequenceSizeOfFlowerGrouping)
    makeTargets(cactusDiskDatabaseString=cactusDiskDatabaseString, configNode=configNode, 
                flowersAndSizes=childFlowers, parentTarget=target, 
                target=childTarget)

class CactusCafDown(CactusRecursionTarget):
    """This target does the get flowers down pass for the CAF alignment phase.
    """    
    def run(self):
        makeChildTargets(self.cactusDiskDatabaseString, self.configNode, self.flowerNames, self, CactusCafDown)
        minSequenceSizeOfFlower = getOptionalAttrib(self.configNode, "minFlowerSize", int, 1)
        maxSequenceSizeOfFlowerGrouping=getOptionalAttrib(self.configNode, "maxFlowerGroupSize", int, default=MAX_SEQUENCE_SIZE_OF_FLOWER_GROUPING)
        childFlowers = runCactusExtendFlowers(self.cactusDiskDatabaseString, self.flowerNames, 
                                              minSequenceSizeOfFlower=minSequenceSizeOfFlower, 
                                              maxSequenceSizeOfFlowerGrouping=maxSequenceSizeOfFlowerGrouping)
        makeTargets(self.cactusDiskDatabaseString, self.configNode, childFlowers, 
                    parentTarget=self, target=CactusCoreWrapper1,
                    overlargeTarget=CactusBlastWrapper)
       
class CactusBlastWrapper(CactusRecursionTarget):
    """Runs blast on the given flower and passes the resulting alignment to cactus core.
    """
    def run(self):
        logger.info("Starting the cactus aligner target")
        #Generate a temporary file to hold the alignments
        alignmentFile = getTempFile(".fa", self.getGlobalTempDir())
        logger.info("Got an alignments file")
        
        #Now make the child aligner target
        blastNode = self.configNode.find("blast")
        blastOptions = \
        makeBlastFromOptions(MakeBlastOptions(chunkSize=getOptionalAttrib(blastNode, "chunkSize", int),
                                              overlapSize=getOptionalAttrib(blastNode, "overlapSize", int),
                                              lastzArguments=getOptionalAttrib(blastNode, "lastzArguments"),
                                              chunksPerJob=getOptionalAttrib(blastNode, "chunksPerJob", int),
                                              compressFiles=getOptionalAttrib(blastNode, "compressFiles", bool)))
        assert len(self.flowerNames) == 2
        assert self.flowerNames[1] == 1
        self.addChildTarget(MakeSequences(self.cactusDiskDatabaseString, 
                                          self.flowerNames[0], alignmentFile, blastOptions=blastOptions,
                                          minimumSequenceLength=getOptionalAttrib(blastNode, "minimumSequenceLengthForBlast", int, 1)))
        logger.info("Created the cactus_aligner child target")
        
        #Now setup a call to cactus core wrapper as a follow on
        self.setFollowOnTarget(CactusCoreWrapper2(self.cactusDiskDatabaseString, self.configNode, [ self.flowerNames, alignmentFile ]))
        logger.info("Setup the follow on cactus_core target")

def runCactusCoreInWorkflow(self, flowerNames, alignmentFile):
    blastParameters = self.configNode.find("blast")
    coreParameters = self.configNode.find("core")
    messages = runCactusCore(cactusDiskDatabaseString=self.cactusDiskDatabaseString,
                      alignments=alignmentFile, 
                      flowerNames=flowerNames,
                      annealingRounds=getOptionalAttrib(coreParameters, "annealingRounds"),  
                      deannealingRounds=getOptionalAttrib(coreParameters, "deannealingRounds"),
                      alignRepeatsAtRound=getOptionalAttrib(coreParameters, "alignRepeatsAtRound", float), 
                      trim=getOptionalAttrib(coreParameters, "trim"),
                      minimumTreeCoverage=getOptionalAttrib(coreParameters, "minimumTreeCoverage", float),
                      blockTrim=getOptionalAttrib(coreParameters, "blockTrim", float),
                      minimumBlockDegree=getOptionalAttrib(coreParameters, "minimumBlockDegree", int), 
                      requiredIngroupFraction=getOptionalAttrib(coreParameters, "requiredIngroupFraction", float),
                      requiredOutgroupFraction=getOptionalAttrib(coreParameters, "requiredOutgroupFraction", float),
                      requiredAllFraction=getOptionalAttrib(coreParameters, "requiredAllFraction", float),
                      singleCopyIngroup=getOptionalAttrib(coreParameters, "singleCopyIngroup", bool),
                      singleCopyOutgroup=getOptionalAttrib(coreParameters, "singleCopyOutgroup", bool),
                      lastzArguments=getOptionalAttrib(blastParameters, "lastzArguments"),
                      minimumSequenceLengthForBlast=getOptionalAttrib(blastParameters, "minimumSequenceLengthForBlast", int, 1),
                      maxAdjacencyComponentSizeRatio=getOptionalAttrib(coreParameters, "maxAdjacencyComponentSizeRatio", float))
    for message in messages:
        self.logToMaster(message)

class CactusCoreWrapper1(CactusRecursionTarget):
    """Runs cactus_core upon a set of flowers and no alignment file.
    """
    def run(self):
        runCactusCoreInWorkflow(self, self.flowerNames, None)
        
class CactusCoreWrapper2(CactusRecursionTarget):
    """Runs cactus_core upon a one flower and one alignment file.
    """
    def run(self):
        runCactusCoreInWorkflow(self, self.flowerNames[0], self.flowerNames[1])
        
############################################################
############################################################
############################################################
#The BAR phase.
#
#Creates the reconstruction structure with blocks
############################################################
############################################################
############################################################

class CactusBarDown(CactusRecursionTarget):
    """This target does the get flowers down pass for the BAR alignment phase.
    """
    def run(self):
        makeChildTargets(self.cactusDiskDatabaseString, self.configNode, self.flowerNames, self, CactusBarDown)
        baseNode = self.configNode.find("base")
        maxSequenceSizeOfFlowerGrouping=getOptionalAttrib(baseNode, "maxFlowerGroupSize", int, default=MAX_SEQUENCE_SIZE_OF_FLOWER_GROUPING)
        childFlowers = runCactusExtendFlowers(self.cactusDiskDatabaseString, self.flowerNames, 
                                              minSequenceSizeOfFlower=1, 
                                              maxSequenceSizeOfFlowerGrouping=maxSequenceSizeOfFlowerGrouping)
        makeTargets(self.cactusDiskDatabaseString, baseNode, childFlowers, parentTarget=self, target=CactusBaseLevelAlignerWrapper, 
                    #cpu=getOptionalAttrib(self.configNode, "numThreads", int, default=sys.maxint),
                    overlargeCpu=getOptionalAttrib(self.configNode, "numThreads", int, default=sys.maxint))

class CactusBaseLevelAlignerWrapper(CactusRecursionTarget):
    """Runs cactus_baseAligner (the BAR algorithm implementation.
    """
    def run(self):
        cpu = 1
        if self.getCpu() != sys.maxint:
            cpu = sys.maxint
        runCactusBaseAligner(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                             flowerNames=self.flowerNames, 
                             maximumLength=getOptionalAttrib(self.configNode, "bandingLimit", float),
                             spanningTrees=getOptionalAttrib(self.configNode, "spanningTrees", int), 
                             gapGamma=getOptionalAttrib(self.configNode, "gapGamma", float), 
                             useBanding=getOptionalAttrib(self.configNode, "useBanding", bool),
                             maxBandingSize=getOptionalAttrib(self.configNode, "maxBandingSize", int),
                             minBandingSize=getOptionalAttrib(self.configNode, "minBandingSize", int), 
                             minBandingConstraintDistance=getOptionalAttrib(self.configNode, "minBandingConstraintDistance", int), 
                             minTraceBackDiag=getOptionalAttrib(self.configNode, "minTraceBackDiag", int), 
                             minTraceGapDiags=getOptionalAttrib(self.configNode, "minTraceGapDiags", int), 
                             constraintDiagonalTrim=getOptionalAttrib(self.configNode, "constaintDiagonalTrim", int), 
                             minimumBlockDegree=getOptionalAttrib(self.configNode, "minimumBlockDegree", int),
                             alignAmbiguityCharacters=getOptionalAttrib(self.configNode, "alignAmbiguityCharacters", bool),
                             pruneOutStubAlignments=getOptionalAttrib(self.configNode, "pruneOutStrubAlignments", bool),
                             requiredIngroupFraction=getOptionalAttrib(self.configNode, "requiredIngroupFraction", float),
                             requiredOutgroupFraction=getOptionalAttrib(self.configNode, "requiredOutgroupFraction", float),
                             requiredAllFraction=getOptionalAttrib(self.configNode, "requiredAllFraction", float),
                             numThreads=cpu)
        
############################################################
############################################################
############################################################
#Normalisation pass
############################################################
############################################################
############################################################
    
class CactusNormalPhase(CactusPhasesTarget):
    def run(self):
        self.logToMaster("Starting the normalisation phase at %s seconds" % time.time())
        normalisationNode = self.options.config.find("normal")
        assert normalisationNode != None
        normalisationIterations = getOptionalAttrib(normalisationNode, "iterations", int, default=1)
        if self.iteration < normalisationIterations:
            self.addChildTarget(CactusNormalDown(self.options.cactusDiskDatabaseString, extractNode(normalisationNode), [ self.flowerName, 1 ]))
        if self.iteration-1 > 0:
            self.setFollowOnTarget(CactusNormalPhase(self.options, self.flowerName, self.iteration+1))
        else:
            self.setFollowOnTarget(CactusPhylogenyPhase(self.options, self.flowerName))
     
class CactusNormalDown(CactusRecursionTarget):
    """This target does the down pass for the normal phase.
    """
    def run(self):
        makeChildTargets(self.cactusDiskDatabaseString, self.configNode, self.flowerNames, self, CactusNormalDown)
        self.setFollowOnTarget(CactusNormalRunnable(self.cactusDiskDatabaseString, self.configNode, self.flowerNames))
        
class CactusNormalRunnable(CactusRecursionTarget):
    """This targets run the normalisation script.
    """ 
    def run(self):
        runCactusMakeNormal(self.cactusDiskDatabaseString, flowerNames=self.flowerNames, 
                            maxNumberOfChains=getOptionalAttrib(self.configNode, "maxNumberOfChains", int))

############################################################
############################################################
############################################################
#Phylogeny pass
############################################################
############################################################
############################################################
    
class CactusPhylogenyPhase(CactusPhasesTarget):        
    def run(self):
        self.logToMaster("Starting the phylogeny phase at %s seconds" % time.time())
        phylogenyNode = self.options.config.find("phylogeny")
        buildTrees = getOptionalAttrib(phylogenyNode, "buildTrees", bool, default=False) or self.options.buildTrees
        if buildTrees:
            self.addChildTarget(CactusPhylogeny(self.options.cactusDiskDatabaseString, extractNode(phylogenyNode), [ self.flowerName, 1 ]))
        self.setFollowOnTarget(CactusReferencePhase(self.options, self.flowerName))

class CactusPhylogeny(CactusRecursionTarget):
    """This target does the down pass for the phylogeny phase.
    """
    def run(self):
        runCactusPhylogeny(self.cactusDiskDatabaseString, flowerNames=self.flowerNames)
        makeChildTargets(self.cactusDiskDatabaseString, self.configNode, self.flowerNames, self, CactusPhylogeny)

############################################################
############################################################
############################################################
#Reference pass
############################################################
############################################################
############################################################
    
class CactusReferencePhase(CactusPhasesTarget):     
    def run(self):
        self.logToMaster("Starting the reference phase at %s seconds" % time.time())
        referenceNode = self.options.config.find("reference")
        buildReference = getOptionalAttrib(referenceNode, "buildReference", bool, default=False) or self.options.buildReference
        if buildReference:
            self.addChildTarget(CactusReferenceDown(self.options.cactusDiskDatabaseString, extractNode(referenceNode), 
                                                    [ self.flowerName, 1 ]))
            self.setFollowOnTarget(CactusReferenceSetCoordinatesUpPhase(self.options, self.flowerName))
        else:
            self.setFollowOnTarget(CactusFacesPhase(self.options, self.flowerName))
        
class CactusReferenceDown(CactusRecursionTarget):
    """This target does the down pass for the reference phase.
    """
    def run(self):
        runCactusReference(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                           flowerNames=self.flowerNames, 
                           matchingAlgorithm=getOptionalAttrib(self.configNode, "matchingAlgorithm"), 
                           permutations=getOptionalAttrib(self.configNode, "permutations", int),
                           referenceEventString=getOptionalAttrib(self.configNode, "reference"), 
                           useSimulatedAnnealing=getOptionalAttrib(self.configNode, "useSimulatedAnnealing", bool),
                           theta=getOptionalAttrib(self.configNode, "theta", float),
                           maxNumberOfChainsBeforeSwitchingToFast=getOptionalAttrib(self.configNode, "maxNumberOfChainsBeforeSwitchingToFast", int))
        makeChildTargets(self.cactusDiskDatabaseString, self.configNode, self.flowerNames, self, CactusReferenceDown)

############################################################
############################################################
############################################################
#Reference coordinates pass, in which coordinates and bases are added to the reference thread
############################################################
############################################################
############################################################

class CactusReferenceSetCoordinatesUpPhase(CactusPhasesTarget):
    def run(self):
        self.logToMaster("Starting the reference coordinate up phase at %s seconds" % time.time())
        referenceNode = self.options.config.find("reference")
        assert referenceNode != None
        self.addChildTarget(CactusSetReferenceCoordinatesUp(self.options.cactusDiskDatabaseString, extractNode(referenceNode), [ self.flowerName, 1 ]))
        self.setFollowOnTarget(CactusSetReferenceCoordinatesDownPhase(self.options, self.flowerName))

class CactusSetReferenceCoordinatesUp(CactusRecursionTarget):
    """Does the up pass for filling Fills in the coordinates, once a reference is added.
    """ 
    def run(self):
        makeChildTargets(self.cactusDiskDatabaseString, self.configNode, self.flowerNames, self, CactusSetReferenceCoordinatesUp)
        self.setFollowOnTarget(CactusSetReferenceCoordinatesUpRunnable(self.cactusDiskDatabaseString, self.configNode, self.flowerNames))
        
class CactusSetReferenceCoordinatesUpRunnable(CactusRecursionTarget):
    """Does the up pass for filling Fills in the coordinates, once a reference is added.
    """ 
    def run(self):
        runCactusAddReferenceCoordinates(cactusDiskDatabaseString=self.cactusDiskDatabaseString, flowerNames=self.flowerNames,
                                         referenceEventString=getOptionalAttrib(self.configNode, "reference"), bottomUpPhase=True)
        
class CactusSetReferenceCoordinatesDownPhase(CactusPhasesTarget):
    def run(self):
        self.logToMaster("Starting the reference coordinate down phase at %s seconds" % time.time())
        referenceNode = self.options.config.find("reference")
        assert referenceNode != None
        self.addChildTarget(CactusSetReferenceCoordinatesDown(self.options.cactusDiskDatabaseString, extractNode(referenceNode), 
                                                              [ self.flowerName, 1 ]))
        self.setFollowOnTarget(CactusFacesPhase(self.options, self.flowerName))
        
class CactusSetReferenceCoordinatesDown(CactusRecursionTarget):
    """Does the down pass for filling Fills in the coordinates, once a reference is added.
    """        
    def run(self):
        makeChildTargets(self.cactusDiskDatabaseString, self.configNode, self.flowerNames, self, CactusSetReferenceCoordinatesDown)
        runCactusAddReferenceCoordinates(cactusDiskDatabaseString=self.cactusDiskDatabaseString, flowerNames=self.flowerNames,
                                         referenceEventString=getOptionalAttrib(self.configNode, "reference"),
                                         bottomUpPhase=False)
        
############################################################
############################################################
############################################################
#Faces pass
############################################################
############################################################ 
############################################################
    
class CactusFacesPhase(CactusPhasesTarget):
    def run(self):
        logger.info("Starting the faces phase")
        facesNode = self.options.config.find("faces")
        buildFaces = getOptionalAttrib(facesNode, "buildFaces", bool, default=False) or self.options.buildFaces
        if buildFaces:
            self.addChildTarget(CactusFaces(self.options.cactusDatabaseString, extractNode(facesNodes), 
                                            [ self.flowerName, 1 ]))
        self.setFollowOnTarget(CactusCheckPhase(self.options, self.flowerName))
        
class CactusFaces(CactusRecursionTarget):
    """This target does the down pass for the faces phase.
    """    
    def run(self):
        runCactusAdjacencies(cactusDiskDatabaseString=self.cactusDiskDatabaseString, flowerNames=self.flowerNames)
        makeChildTargets(self.cactusDiskDatabaseString, self.configNode, self.flowerNames, self, CactusFaces)

############################################################
############################################################
############################################################
#Check pass
############################################################
############################################################
############################################################
    
class CactusCheckPhase(CactusPhasesTarget):
    def run(self):
        checkNode=self.options.config.find("check")
        assert checkNode != None
        if getOptionalAttrib(checkNode, "runCheck", bool, default=False): #self.options.skipCheck and 0:
            self.logToMaster("Starting the verification phase at %s seconds" % (time.time()))
            self.addChildTarget(CactusCheck(self.options.cactusDiskDatabaseString, extractNode(checkNode), 
                                            [ self.flowerName, 1 ]))
        self.setFollowOnTarget(CactusMafGeneratorPhase(self.options, self.flowerName))
        
class CactusCheck(CactusRecursionTarget):
    """This target does the down pass for the check phase.
    """
    def run(self):
        runCactusCheck(self.cactusDiskDatabaseString, self.flowerNames, getOptionalAttrib(self.configNode, "checkNormalised", bool))
        makeChildTargets(self.cactusDiskDatabaseString, self.configNode, self.flowerNames, self, CactusCheck)   

############################################################
############################################################
############################################################
#Check pass
############################################################
############################################################
############################################################

class CactusMafGeneratorPhase(CactusPhasesTarget):
    def run(self):
        self.logToMaster("Starting the maf generation phase at %s seconds" % time.time())
        if self.options.buildReference: #Must have reference building set.
            referenceNode = self.options.config.find("reference")
            mafGeneratorNode = self.options.config.find("maf")
            if referenceNode.attrib.has_key("reference"):
                mafGeneratorNode.attrib["reference"] = referenceNode.attrib["reference"]
            if getOptionalAttrib(mafGeneratorNode, "buildMaf", bool, default=False) or self.options.buildMaf:
                self.addChildTarget(CactusMafGeneratorUp(cactusDiskDatabaseString=self.options.cactusDiskDatabaseString, 
                                                         configNode=extractNode(mafGeneratorNode), 
                                                         flowerNames=[ self.flowerName, 1 ], 
                                                         parentTempDir=None, outputFile=self.options.experimentFile.find("maf").attrib["path"]))

class CactusMafGeneratorUp(CactusRecursionTarget):
    """Generate the maf my merging mafs from the children.
    """ 
    def __init__(self, cactusDiskDatabaseString, configNode, flowerNames, parentTempDir, outputFile, memory=sys.maxint, cpu=sys.maxint):
        CactusRecursionTarget.__init__(self, cactusDiskDatabaseString, configNode, flowerNames, memory, cpu)
        self.parentTempDir = parentTempDir
        self.outputFile = outputFile
    
    def run(self):
        def fn(cactusDiskDatabaseString, configNode, flowerNames, memory=sys.maxint, cpu=sys.maxint):
            return CactusMafGeneratorUp(cactusDiskDatabaseString=cactusDiskDatabaseString, configNode=configNode, 
                                        flowerNames=flowerNames, parentTempDir=self.getGlobalTempDir(), 
                                        outputFile=None,
                                        memory=memory, cpu=cpu)
        #(self, cactusDiskDatabaseString, configNode, flowerNames, memory=sys.maxint, cpu=sys.maxint):
        makeChildTargets(self.cactusDiskDatabaseString, self.configNode, self.flowerNames, self, fn)
        self.setFollowOnTarget(CactusMafGeneratorUpRunnable(self.cactusDiskDatabaseString, 
                                                            self.configNode, self.flowerNames, self.parentTempDir, self.outputFile))
        
class CactusMafGeneratorUpRunnable(CactusMafGeneratorUp):
    """Does the up pass for filling in the coordinates, once a reference is added.
    """ 
    def run(self):
        runCactusRecursiveMafGenerator(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                              flowerNames=self.flowerNames,
                              referenceEventString=getOptionalAttrib(self.configNode, "reference"), #self.configNode.attrib["reference"], #getOptionalAttrib(self.configNode, "reference"), 
                              childDir=self.getGlobalTempDir(), parentDir=self.parentTempDir, 
                              outputFile=self.outputFile,
                              showOnlySubstitutionsWithRespectToReference=\
                              getOptionalAttrib(self.configNode, "showOnlySubstitutionsWithRespectToReference", bool))

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
    
    parser.add_option("--buildMaf", dest="buildMaf", action="store_true",
                      help="Build a maf", default=False)
    
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
