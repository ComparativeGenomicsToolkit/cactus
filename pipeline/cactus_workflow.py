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

from cactus.blastAlignment.cactus_aligner import MakeSequences
from cactus.blastAlignment.cactus_batch import MakeBlastOptions
from cactus.blastAlignment.cactus_batch import makeBlastFromOptions

from cactus.preprocessor.cactus_preprocessor import BatchPreprocessor

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

def getOptionalAttrib(node, attribName, type=None, default=None):
    """Get an optional attrib, or None, if not set.
    """
    if node.attrib.has_key(attribName):
        if type != None:
            if type == bool:
                type(int(node.attrib[attribName]))
            return type(node.attrib[attribName])
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
        if int(alignmentNode.find("blast_misc").attrib["filterByIdentity"]):
            longestPath = getLongestPath(newickTreeParser(self.options.speciesTree))
            adjustedPath = float(alignmentNode.find("blast_misc").attrib["identityRatio"]) * longestPath + float(alignmentNode.find("blast_misc").attrib["minimumDistance"])
            identity = str(100 - int(100 * inverseJukesCantor(adjustedPath)))
            logger.info("The blast stage will filter by identity, the calculated minimum identity is %s from a longest path of %s and an adjusted path of %s" % (identity, longestPath, adjustedPath))
            for iterationNode in alignmentNode.find("iterations").findall("iteration"):
                if iterationNode.attrib["type"] == "blast":
                    blastNode = iterationNode.find("blast")
                    assert "IDENTITY" in blastNode.attrib["blastString"]
                    blastNode.attrib["blastString"] = blastNode.attrib["blastString"].replace("IDENTITY", identity)
                    assert "IDENTITY" in blastNode.attrib["selfBlastString"]
                    blastNode.attrib["selfBlastString"] = blastNode.attrib["selfBlastString"].replace("IDENTITY", identity)
                else:
                    assert iterationNode.attrib["type"] == "base"

    def run(self):
        self.logToMaster("Starting setup phase target at %s seconds" % time.time())
        #Modify the config options
        self.modifyConfig()
        #Make the child setup job.
        self.addChildTarget(CactusPreprocessorPhase(self.options, self.sequences))
        #initialise the down pass as the follow on.. using special '0'
        self.setFollowOnTarget(CactusAlignmentPhase('0', self.options))
        logger.info("Created child target preprocessor job, and follow on down pass job")

class CactusPreprocessorPhase(Target):
    def __init__(self, options, sequences):
        Target.__init__(self, time=0.0002)
        self.options = options 
        self.sequences = sequences

    def run(self):
        self.logToMaster("Starting preprocessor phase target at %s seconds" % time.time())
        
        processedSequences = self.sequences
       
        prepNode = self.options.config.find("preprocessor")
        if prepNode is not None:
            tempDir = getTempDirectory(self.getGlobalTempDir())
            processedSequences = map(lambda x: tempDir + "/" + x, self.sequences)
            logger.info("Adding child batch_preprocessor target")
            self.addChildTarget(BatchPreprocessor(self.options, self.sequences, self.sequences, tempDir, 0))

        self.setFollowOnTarget(CactusSetupWrapper(self.options, processedSequences))
        logger.info("Created followOn target cactus_setup job, and follow on down pass job")
        
class CactusSetupWrapper(Target):
    def __init__(self, options, sequences):
        Target.__init__(self, time=1.0)
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

def compressFlowerList(flowerList):
    """Make a compressed list of flowers
    """
    return bz2.compress(" ".join([ str(i) for i, j in flowerList ]))

def decompressFlowerList(flowerList):
    """Decompress the list of flowers
    """
    return [ int(i) for i in bz2.decompress(flowerList).split() ]

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
            iterationNode = extractNode(iterations[self.iteration])
            if iterationNode.attrib["type"] == "blast":
                self.addChildTarget(CactusCafDown(self.cactusDiskDatabaseString, iterationNode, compressFlowersList([ self.flowerName, 1 ])))
            else:
                assert iterationNode.attrib["type"] == "base"
                self.addChildTarget(CactusBarDown(self.cactusDiskDatabaseString, iterationNode, compressFlowersList([ self.flowerName, 1 ])))
            self.setFollowOnTarget(CactusAlignmentPhase(self.options, self.flowerName, self.iteration+1))
        else:
            self.setFollowOnTarget(CactusNormalPhase(self.options, self.flowerName, ))

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
    def __init__(self, cactusDiskDatabaseString, configNode, flowerNames):
        targetStringId = str(self.__class__)[8:-2]
        Target.__init__(self, time=getOptionalAttrib(configNode, "time_" + targetStringId, float), 
                        cpu=getOptionalAttrib(configNode, "cpu_" + targetStringId, int), 
                        memory=getOptionalAttrib(configNode, "memory_" + targetStringId, int))
        self.cactusDiskDatabaseString = cactusDiskDatabaseString
        self.configNode = configNode
        self.flowerNames = flowerNames
  
MAX_SEQUENCE_SIZE_OF_FLOWER_GROUPING=10000        

def makeTargets(cactusDiskDatabaseString, configNode, flowersAndSizes, 
                parentTarget, target):
    """Make a set of targets for a given set of flowers.
    """
    flowerNames = []
    totalSequenceSize = 0.0
    maxSequenceSizeOfFlowerGrouping=getOptionalAttrib(configNode, "maxFlowerGroupSize", int, default=MAX_SEQUENCE_SIZE_OF_FLOWER_GROUPING)
    
    for totalFlowerSize, firstFlowerName, flowerNumber in flowersAndSizes:
        if totalFlowerSize > maxSequenceSizeOfFlowerGrouping: #Make sure large flowers are on there own, in their own job
            if flowerNumber != 1:
                logger.critical("Got more than one flower: %s %s %s %s %s" % (firstFlowerName, flowerNumber, totalFlowerSize, len(flowersAndSizes), flowersAndSizes))
            assert flowerNumber == 1
            parentTarget.logToMaster("Adding an oversize flower: %s on its own, with %s bases for target class %s" \
                                     % (firstFlowerName, totalFlowerSize, target))
            parentTarget.addChildTarget(target(cactusDiskDatabaseString=cactusDiskDatabaseString, configNode=configNode, 
                                                       flowerNames=compressFlowersList([ firstFlowerName, 1 ])))
        else:
            totalSequenceSize += totalFlowerSize
            flowerNames.append(firstFlowerName)
            flowerNames.append(flowerNumber)
            if totalSequenceSize >= maxSequenceSizeOfFlowerGrouping: 
                parentTarget.addChildTarget(target(cactusDiskDatabaseString=cactusDiskDatabaseString, configNode=configNode, 
                                                   flowerNames=compressFlowersList(flowerNames)))
                flowerNames = []
                totalSequenceSize = 0.0
    if len(flowerNames) > 0:
        parentTarget.addChildTarget(target(cactusDiskDatabaseString, configNode, compressFlowersList(flowerNames)))
     
def makeChildTargets(cactusDiskDatabaseString, configNode, flowerNames, target, childTarget):
    """Make a set of child targets for a given set of parent flowers.
    """
    minSequenceSizeOfFlower = getOptionalAttrib(configNode, "minFlowerSize", int) 
    maxSequenceSizeOfFlowerGrouping=getOptionalAttrib(configNode, "maxFlowerGroupSize", int, default=MAX_SEQUENCE_SIZE_OF_FLOWER_GROUPING)
    childFlowers = runCactusGetFlowers(options.cactusDiskDatabaseString, flowerNames, 
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
        minSequenceSizeOfFlower = getOptionalAttrib(configNode, "minFlowerSize", int)
        maxSequenceSizeOfFlower=getOptionalAttrib(configNode, "maxFlowerSize", int)  
        maxSequenceSizeOfFlowerGrouping=getOptionalAttrib(configNode, "maxFlowerGroupSize", int, default=MAX_SEQUENCE_SIZE_OF_FLOWER_GROUPING)
        for totalFlowerSize, firstChildName, flowerNumber in runCactusExtendFlowers(self.cactusDiskDatabaseString, 
                                                                                    self.flowerNames, 
                                                                       minSequenceSizeOfFlower=minFlowerSize, 
                                                                       maxSequenceSizeOfFlower=maxFlowerSize,
                                                                       maxSequenceSizeOfFlowerGrouping=maxSequenceSizeOfFlowerGrouping):
            for i in xrange(flowerNumber):
                self.addChildTarget(CactusBlastWrapper(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                                                       configNode=self.configNode, flowerNames=str(firstChildName + i)))
       
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
                                              blastString=getOptionalAttrib(blastNode, "blastString"),
                                              selfBlastString=getOptionalAttrib(blastNode, "selfBlastString"),
                                              chunksPerJob=getOptionalAttrib(blastNode, "chunksPerJob", int),
                                              compressFiles=getOptionalAttrib(blastNode, "compressFiles", bool)))
        self.addChildTarget(MakeSequences(cactusDiskDatabaseString=self.cactusDiskDatabaseString, 
                                          flowerName=self.flowerNames, alignmentFile=alignmentFile, blastOptions=blastOptions,
                                          minimumSequenceLength=getOptionalAttrib(blastNode, "minimumSequenceLength", int, 0)))
        logger.info("Created the cactus_aligner child target")
        
        #Now setup a call to cactus core wrapper as a follow on
        self.setFollowOnTarget(CactusCoreWrapper(self.cactusDiskDatabaseString, self.configNode, [ self.flowerNames, alignmentFile ]))
        logger.info("Setup the follow on cactus_core target")

class CactusCoreWrapper(CactusRecursionTarget):
    """Runs cactus_core upon a given flower and alignment file.
    """
    def run(self):
        coreParameters = self.configNode.find("core")
        runCactusCore(cactusDiskDatabaseString=self.options.cactusDiskDatabaseString,
                      alignmentFile=self.alignmentFile[1], 
                      flowerName=self.flowerNames[0],
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
                      singleCopyOutgroup=getOptionalAttrib(coreParameters, "singleCopyOutgroup", bool))
        
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
        minSequenceSizeOfFlower = getOptionalAttrib(configNode, "minFlowerSize", int)
        maxSequenceSizeOfFlowerGrouping=getOptionalAttrib(configNode, "maxFlowerGroupSize", int, default=MAX_SEQUENCE_SIZE_OF_FLOWER_GROUPING)
        childFlowers = runCactusExtendFlowers(self.cactusDiskDatabaseString, self.flowerNames, 
                                              minSequenceSizeOfFlower=1, 
                                              maxSequenceSizeOfFlowerGrouping=maxSequenceSizeOfFlowerGrouping)
        makeTargets(self.cactusDiskDatabaseString, self.configNode, childFlowers, parentTarget=self, target=CactusBaseLevelAlignerWrapper,
                    maxSequenceSizeOfFlowerGrouping=maxSequenceSizeOfFlowerGrouping)

class CactusBaseLevelAlignerWrapper(CactusRecursionTarget):
    """Runs cactus_baseAligner (the BAR algorithm implementation.
    """
    def run(self):
        assert self.configNode.attrib["type"] == "base"
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
                             numThreads=getOptionalAttrib(self.configNode, "numThreads", int))
        
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
        normalisationNode = self.options.config.find("normalisation")
        assert normalisationNode != None
        normalisationRounds = getOptionalAttrib(normalisationNode, "rounds", int, default=1)
        if self.iteration < self.normalisationRounds:
            self.addChildTarget(CactusNormalDown(self.options.cactusDiskDatabaseString, normalisationNode, compressFlowersList([ self.flowerName, 1 ])))
        if self.normalisationRounds-1 > 0:
            self.setFollowOnTarget(CactusNormalPhase(self.options, self.flowerName, self.iteration+1))
        else:
            self.setFollowOnTarget(CactusPhylogenyPhase(self.options, self.flowerName))
     
class CactusNormalDown(CactusRecursionTarget):
    """This target does the down pass for the normal phase.
    """
    def run(self):
        makeChildTargets(self.cactusDiskDatabaseString, self.configNode, self.flowerNames, self, CactusNormalDown)
        self.setFollowOnTarget(CactusNormalRunnable(options=self.options, flowerNames=self.flowerNames))
        
class CactusNormalRunnable(CactusRecursionTarget):
    """This targets run the normalisation script.
    """ 
    def run(self):
        runCactusMakeNormal(self.cactusDiskDatabaseString, flowerNames=self.flowerNames, 
                            getOptionalAttrib(self.configNode, "maxNumberOfChains", int))

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
        buildTrees = getOptionalAttrib(buildTrees, "buildTrees", bool, default=False)
        if buildTrees:
            self.addChildTarget(CactusPhylogeny(self.options.cactusDiskDatabaseString, phylogenyNode, compressFlowersList([ self.flowerName, 1 ])))
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
        buildReference = getOptionalAttrib(referenceNode, "buildReference", bool, default=False)
        if buildReference:
            self.addChildTarget(CactusReferenceDown(self.options.cactusDiskDatabaseString, referenceNode, 
                                                    compressFlowersList([ self.flowerName, 1 ])))
            self.setFollowOnTarget(CactusReferenceSetCoordinatesUpPhase(self.options, self.flowerName))
        else:
            self.setFollowOnTarget(CactusFacesPhase(self.options, self.flowerName))
        
class CactusReferenceDown(CactusRecursionTarget):
    """This target does the down pass for the reference phase.
    """
    def run(self):
        runCactusReference(cactusDiskDatabaseString=self.options.cactusDiskDatabaseString, 
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
        self.addChildTarget(CactusSetReferenceCoordinatesUp(self.options, None, [ self.flowerName, 1 ]))
        self.setFollowOnTarget(CactusSetReferenceCoordinatesDownPhase(self.flowerName, self.options))

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
        self.addChildTarget(CactusSetReferenceCoordinatesDown(self.cactusDiskDatabaseString, self.configNode, compressFlowersList([ self.flowerName, 1 ])))
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
        buildFaces = getOptionalAttrib(referenceNode, "buildFaces", bool, default=False)
        if buildFaces:
            self.addChildTarget(CactusFaces(self.options.cactusDatabaseString, facesNodes, compressFlowersList([ self.flowerName, 1 ])))
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
        if getOptionalAttrib(checkNode, "runCheck", bool, default=False): #self.options.skipCheck and 0:
            self.logToMaster("Starting the verification phase at %s seconds" % (time.time()))
            self.addChildTarget(CactusCheck(self.options.cactusDiskDatabaseString, checkNode, compressFlowersList([ self.flowerName, 1 ])))
        
class CactusCheck(CactusRecursionTarget):
    """This target does the down pass for the check phase.
    """
    def run(self):
        runCactusCheck(self.options.cactusDiskDatabaseString, self.flowerNames, getOptionalAttrib(self.configNode, "checkNormalised", bool))
        makeChildTargets(self.cactusDiskDatabaseString, self.configNode, self.flowerNames, self, CactusCheck)   

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
    
    options, args = parser.parse_args()
    setLoggingFromOptions(options)
    
    if len(args) != 0:
        raise RuntimeError("Unrecognised input arguments: %s" % " ".join(args))

    # process the options
    expandWorkflowOptions(options)
    #Get the sequences
    sequences = options.experimentFile.attrib["sequences"].split()

    baseTarget = CactusSetupPhase(options, sequences)
    
    Stack(baseTarget).startJobTree(options)
    logger.info("Done with job tree")

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    from cactus.pipeline.cactus_workflow import *
    _test()
    main()
